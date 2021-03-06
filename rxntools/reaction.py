#!/usr/bin/env python
# -*- coding: utf-8 -*-
import copy
from .template import *
from itertools import permutations


class ChemicalReaction(Base):
    def __init__(self, reaction_smarts, special_groups=[]):
        super().__init__(reaction_smarts=reaction_smarts, is_template=False)
        self.Sanitized()
        self.RearangeMols()
        self.GetChiralityInfo()
        self.special_groups = special_groups

    # Functions used for extract reaction template.
    def ExtractTemplate(self, validate=True):
        changed_atom_tags = list(map(str, self.reacting_atoms_mn))
        reactant_fragments = self.__GetFragmentsForReactingAtoms(
            category='reactant', depth=1
        )
        expansion = self.__ExpandChangedAtomTags(reactant_fragments)
        product_fragments = self.__GetFragmentsForReactingAtoms(
            category='product', depth=0, expansion=expansion
        )
        template_smarts = '{}>>{}'.format(reactant_fragments, product_fragments)
        template = ReactionTemplate(template_smarts)
        if validate:
            if template.rxn.Validate()[1] != 0:
                raise RuntimeError(
                    f'Could not validate reaction successfully. \n'
                    f'Input reaction smarts: {self.reaction_smarts}.\n'
                    f'Extracted reaction templates: {template_smarts}')
            reactants = copy.copy(self.reactants)
            self._RemoveAtomMap(reactants)
            SMILES_r = list(map(Chem.MolToSmiles, reactants))
            reactants = list(map(Chem.MolFromSmiles, SMILES_r))
            products = copy.copy(self.products)
            self._RemoveAtomMap(products)
            SMILES_p = Chem.MolToSmiles(CombineMols(products))
            SMILES_template_p = []
            for reactants_ in list(permutations(reactants, len(reactants))):
                template_products = template.rxn.RunReactants(reactants_)
                SMILES_template_p += [Chem.MolToSmiles(CombineMols(products))
                                      for products in template_products]
            if not SMILES_p in SMILES_template_p:
                raise RuntimeError(
                    f'For chemical reaction:\n{self.reaction_smarts}\n'
                    f'and extracted reaction template:\n{template.Smarts}\n'
                    f'The true products:\n{SMILES_p}\n'
                    f'not in the products obtained by runing template\n{SMILES_template_p}'
                )
            assert (SMILES_p in SMILES_template_p)
        return template.Smarts

    def __GetFragmentsForReactingAtoms(self, category, depth, expansion=[]):
        mols = self.reactants if category == 'reactant' else self.products
        fragments = []
        for mol in mols:
            symbol_replacements = []

            if category == 'reactant':
                groups = self.__GetSpecialGroups(mol)
            else:
                groups = []

            atoms_to_use = []
            for atom in mol.GetAtoms():
                # Check self (only tagged atoms)
                if getAtomMapNumber(atom) in self.ReactingAtomsMN:
                    atoms_to_use.append(atom.GetIdx())
                    symbol = self.__GetAtomSmarts(atom)
                    if symbol != atom.GetSmarts():
                        symbol_replacements.append((atom.GetIdx(), symbol))
                    continue

            # Check neighbors (any atom) and special groups in reactants
            # Inactive when set radius=0
            for k in range(depth):
                atoms_to_use, symbol_replacements = self.__ExpandAtomsToSse(
                    mol,
                    atoms_to_use,
                    groups=groups,
                    symbol_replacements=symbol_replacements
                )

            if category == 'product':
                # Add extra labels to include (for products only)
                if expansion:
                    for atom in mol.GetAtoms():
                        if ':' not in atom.GetSmarts():
                            continue
                        label = atom.GetSmarts().split(':')[1][:-1]
                        if label in expansion and label not in self.ReactingAtomsMN:
                            atoms_to_use.append(atom.GetIdx())
                            # Make the expansion a wildcard
                            symbol_replacements.append(
                                (atom.GetIdx(), self.__ConvertAtomToWildcard(atom)))

            # Define new symbols to replace terminal species with wildcards
            # (don't want to restrict templates too strictly)
            symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
            for (i, symbol) in symbol_replacements:
                symbols[i] = symbol

            if not atoms_to_use:
                continue

            fragments.append('(' + Chem.MolFragmentToSmiles(
                mol, atoms_to_use,
                atomSymbols=symbols,
                allHsExplicit=True,
                isomericSmiles=True,
                allBondsExplicit=True) + ')')
        fragments.sort()
        return '.'.join(fragments)

    def __ExpandAtomsToSse(self, mol, atoms_to_use, groups,
                           symbol_replacements):
        '''Copy from Rdchiral'''

        # Copy
        new_atoms_to_use = atoms_to_use[:]
        # Look for all atoms in the current list of atoms to use
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in atoms_to_use:
                continue
            # Ensure membership of changed atom is checked against group
            for group in groups:
                if int(atom.GetIdx()) in group[0]:
                    for idx in group[1]:
                        if idx not in atoms_to_use:
                            new_atoms_to_use.append(idx)
                            symbol_replacements.append(
                                (idx, self.__ConvertAtomToWildcard(
                                    mol.GetAtomWithIdx(idx))))
            # Look for all nearest neighbors of the currently-included atoms
            for neighbor in atom.GetNeighbors():
                # Evaluate nearest neighbor atom to determine what should be included
                new_atoms_to_use, symbol_replacements = \
                    self.__ExpandAtomsToUseAtom(
                        new_atoms_to_use,
                        neighbor,
                        groups=groups,
                        symbol_replacements=symbol_replacements)

        return new_atoms_to_use, symbol_replacements

    def __ExpandAtomsToUseAtom(self, atoms_to_use, atom, groups,
                               symbol_replacements):
        '''Copy from Rdchiral'''
        found_in_group = False
        for group in groups:  # first index is atom IDs for match, second is what to include
            if atom.GetIdx() in group[0]:  # int correction
                # Add the whole list, redundancies don't matter
                # *but* still call convert_atom_to_wildcard!
                for idx in group[1]:
                    if idx not in atoms_to_use:
                        atoms_to_use.append(idx)
                        symbol_replacements.append(
                            (idx, self.__ConvertAtomToWildcard(atom)))
                found_in_group = True
        if found_in_group:
            return atoms_to_use, symbol_replacements

        # How do we add an atom that wasn't in an identified important functional group?
        # Develop generalized SMARTS symbol

        # Skip current candidate atom if it is already included
        if atom.GetIdx() in atoms_to_use:
            return atoms_to_use, symbol_replacements

        # Include this atom
        atoms_to_use.append(atom.GetIdx())

        # Look for suitable SMARTS replacement
        symbol_replacements.append(
            (atom.GetIdx(),
             self.__ConvertAtomToWildcard(atom)))

        return atoms_to_use, symbol_replacements

    def __ExpandChangedAtomTags(self, reactant_fragments):
        '''Copy from Rdchiral'''

        expansion = []
        atom_tags_in_reactant_fragments = re.findall('\:([0-9]+)\]',
                                                     reactant_fragments)
        for atom_tag in atom_tags_in_reactant_fragments:
            if int(atom_tag) not in self.ReactingAtomsMN:
                expansion.append(atom_tag)
        return expansion

    def __GetSpecialGroups(self, mol):
        '''Copy from Rdchiral'''
        if self.special_groups == 'rdchiral':
            # Define templates
            group_templates = [
                (range(3), '[OH0,SH0]=C[O,Cl,I,Br,F]',),
                # carboxylic acid / halogen
                (range(3), '[OH0,SH0]=CN',),  # amide/sulfamide
                (range(4), 'S(O)(O)[Cl]',),  # sulfonyl chloride
                (range(3), 'B(O)O',),  # boronic acid/ester
                ((0,), '[Si](C)(C)C'),  # trialkyl silane
                ((0,), '[Si](OC)(OC)(OC)'),  # trialkoxy silane, default to methyl
                (range(3), '[N;H0;$(N-[#6]);D2]-,=[N;D2]-,=[N;D1]',),  # azide
                (range(8), 'O=C1N([Br,I,F,Cl])C(=O)CC1',),  # NBS brominating agent
                (range(11), 'Cc1ccc(S(=O)(=O)O)cc1'),  # Tosyl
                ((7,), 'CC(C)(C)OC(=O)[N]'),  # N(boc)
                ((4,), '[CH3][CH0]([CH3])([CH3])O'),  #
                (range(2), '[C,N]=[C,N]',),  # alkene/imine
                (range(2), '[C,N]#[C,N]',),  # alkyne/nitrile
                ((2,), 'C=C-[*]',),  # adj to alkene
                ((2,), 'C#C-[*]',),  # adj to alkyne
                ((2,), 'O=C-[*]',),  # adj to carbonyl
                ((3,), 'O=C([CH3])-[*]'),  # adj to methyl ketone
                ((3,), 'O=C([O,N])-[*]',),  # adj to carboxylic acid/amide/ester
                (range(4), 'ClS(Cl)=O',),  # thionyl chloride
                (range(2), '[Mg,Li,Zn,Sn][Br,Cl,I,F]',),
                # grinard/metal (non-disassociated)
                (range(3), 'S(O)(O)',),  # SO2 group
                (range(2), 'N~N',),  # diazo
                ((1,), '[!#6;R]@[#6;R]',),  # adjacency to heteroatom in ring
                ((2,), '[a!c]:a:a',),
                # two-steps away from heteroatom in aromatic ring
                # ((1,), 'c(-,=[*]):c([Cl,I,Br,F])',), # ortho to halogen on ring - too specific?
                # ((1,), 'c(-,=[*]):c:c([Cl,I,Br,F])',), # meta to halogen on ring - too specific?
                ((0,), '[B,C](F)(F)F'),  # CF3, BF3 should have the F3 included
            ]

            # Stereo-specific ones (where we will need to include neighbors)
            # Tetrahedral centers should already be okay...
            group_templates += [
                ((1, 2,), '[*]/[CH]=[CH]/[*]'),  # trans with two hydrogens
                ((1, 2,), '[*]/[CH]=[CH]\[*]'),  # cis with two hydrogens
                ((1, 2,), '[*]/[CH]=[CH0]([*])\[*]'),  # trans with one hydrogens
                ((1, 2,), '[*]/[D3;H1]=[!D1]'),
                # specified on one end, can be N or C
            ]

            # Build list
            groups = []
            for (add_if_match, template) in group_templates:
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts(template),
                                                  useChirality=True)
                for match in matches:
                    add_if = []
                    for pattern_idx, atom_idx in enumerate(match):
                        if pattern_idx in add_if_match:
                            add_if.append(atom_idx)
                    groups.append((add_if, match))
            return groups
        else:
            return self.special_groups

    @staticmethod
    def __ConvertAtomToWildcard(atom):
        '''copy from Rdchiral'''

        # Is this a terminal atom? We can tell if the degree is one
        if atom.GetDegree() == 1:
            symbol = '[' + atom.GetSymbol() + ';D1;H{}'.format(
                atom.GetTotalNumHs())
            if atom.GetFormalCharge() != 0:
                charges = re.search('([-+]+[1-9]?)', atom.GetSmarts())
                symbol = symbol.replace(';D1', ';{};D1'.format(charges.group()))

        else:
            # Initialize
            symbol = '['

            # Add atom primitive - atomic num and aromaticity (don't use COMPLETE wildcards)
            if atom.GetAtomicNum() != 6:
                symbol += '#{};'.format(atom.GetAtomicNum())
                if atom.GetIsAromatic():
                    symbol += 'a;'
            elif atom.GetIsAromatic():
                symbol += 'c;'
            else:
                symbol += 'C;'

            # Charge is important
            if atom.GetFormalCharge() != 0:
                charges = re.search('([-+]+[1-9]?)', atom.GetSmarts())
                if charges: symbol += charges.group() + ';'

            # Strip extra semicolon
            if symbol[-1] == ';': symbol = symbol[:-1]

        # Close with label or with bracket
        label = re.search('\:[0-9]+\]', atom.GetSmarts())
        if label:
            symbol += label.group()
        else:
            symbol += ']'

        return symbol

    @staticmethod
    def __GetAtomSmarts(atom):
        symbol = atom.GetSmarts()
        # Sometimes the Aromaticity of a atom is changed after sanitize, but its
        # Smarts didnt updated.
        if not atom.GetIsAromatic() and symbol[1].islower():
            symbol = symbol[0] + symbol[1].upper() + symbol[2:]
        # CUSTOM SYMBOL CHANGES
        if atom.GetTotalNumHs() == 0:
            # Be explicit when there are no hydrogens
            if ':' in symbol:  # stick H0 before label
                symbol = symbol.replace(':', ';H0:')
            else:  # stick before end
                symbol = symbol.replace(']', ';H0]')

        # print('Being explicit about H0!!!!')
        if atom.GetFormalCharge() == 0:
            # Also be explicit when there is no charge
            if ':' in symbol:
                symbol = symbol.replace(':', ';+0:')
            else:
                symbol = symbol.replace(']', ';+0]')
        return symbol

    # Functions to automatically correct reactions.
    def FixChiralityInfo(self, rule):
        if rule == 'complete':
            self.__CompleteMissingChirality()
        elif rule == 'delete':
            self.__DeleteMissingChirality()
        else:
            raise RuntimeError(f'Unknown chirality fix rule {rule}')
        return ChemicalReaction(self.Smarts)

    def __CompleteMissingChirality(self):
        """
        If a chiral atom find in from_record but not in to_record, and their
        neighbors are the same, then the chiraly record is copyed.
        :param from_records:
        :param to_records:
        :param category: reactant or product
        :return:
        """
        self.__ModifyMissingChirality(
            self.reactants_chirality_record,
            self.products_chirality_record,
            'product',
            'complete'
        )
        self.__ModifyMissingChirality(
            self.products_chirality_record,
            self.reactants_chirality_record,
            'reactant',
            'complete'
        )

    def __DeleteMissingChirality(self):
        self.__ModifyMissingChirality(
            self.reactants_chirality_record,
            self.products_chirality_record,
            'product',
            'delete'
        )
        self.__ModifyMissingChirality(
            self.products_chirality_record,
            self.reactants_chirality_record,
            'reactant',
            'delete'
        )

    def __ModifyMissingChirality(self, records1, records2, category2, rule):
        mols = self.reactants if category2 == 'reactant' else self.products
        delete_keys = []
        for mn, record in records1.items():
            if mn not in records2:
                atom_unchiral = self._GetAtomFromMapNumber(
                    mols, mn)
                acr = self._GetAtomChiralityRecord(atom_unchiral,
                                                   category2)
                if '.'.join(list(map(str, record[1]))) == \
                        '.'.join(list(map(str, acr[2]))):
                    if rule == 'complete':
                        records2[mn] = record
                    elif rule == 'delete':
                        delete_keys.append(mn)
        for key in delete_keys:
            records1.pop(key)
