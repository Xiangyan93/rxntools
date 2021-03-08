#!/usr/bin/env python
# -*- coding: utf-8 -*-
import warnings
import copy
import numpy as np
from rdkit.Chem import rdChemReactions
from.substructure import *


class Base:
    def __init__(self, reaction_smarts, is_template):
        self.reaction_smarts = reaction_smarts
        self.is_template = is_template
        self.rxn = Chem.ReactionFromSmarts(reaction_smarts)
        self.reacting_atoms_mn = None

    @property
    def reactants(self):
        return list(self.rxn.GetReactants())

    @property
    def reactants_SMILES(self):
        if not self.reactants:
            return None
        reactants = copy.deepcopy(self.reactants)
        self._RemoveAtomMap(reactants)
        SMILES = Chem.MolToSmiles(CombineMols(reactants))
        # SMILES with '[C]' will be problematic when applying template on it
        SMILES = re.sub('\[([CNOcno])]', lambda x: x.group(1), SMILES)
        return Chem.MolToSmiles(Chem.MolFromSmiles(SMILES))

    @property
    def products(self):
        return list(self.rxn.GetProducts())

    @property
    def products_SMILES(self):
        if not self.products:
            return None
        products = copy.deepcopy(self.products)
        self._RemoveAtomMap(products)
        SMILES = Chem.MolToSmiles(CombineMols(products))
        # SMILES with '[C]' will be problematic when applying template on it
        SMILES = re.sub('\[([CNOcno])]', lambda x: x.group(1), SMILES)
        return Chem.MolToSmiles(Chem.MolFromSmiles(SMILES))

    @property
    def agents(self):
        return list(self.rxn.GetAgents())

    @property
    def Smarts(self):
        if self.is_template:
            reactants_smarts = self._Frag2CanonicalSmarts(self.reactants)
            reactants_smarts = self._AddChiralOnSmarts(
                reactants_smarts, self.reactants_chirality_record)
            products_smarts = self._Frag2CanonicalSmarts(self.products)
            products_smarts = self._AddChiralOnSmarts(
                products_smarts, self.products_chirality_record)
            return reactants_smarts + '>' + \
                   self._Mols2CanonicalSmiles(self.agents) + '>' + \
                   products_smarts
        else:
            reactants_smarts = self._Mols2CanonicalSmarts(self.reactants)
            reactants_smarts = self._AddChiralOnSmarts(
                reactants_smarts, self.reactants_chirality_record)
            products_smarts = self._Mols2CanonicalSmarts(self.products)
            products_smarts = self._AddChiralOnSmarts(
                products_smarts, self.products_chirality_record)
            return reactants_smarts + '>' + \
                   self._Mols2CanonicalSmiles(self.agents) + '>' + \
                   products_smarts

    @property
    def ReactingAtomsMN(self):
        if self.reacting_atoms_mn is None:
            self.reacting_atoms_mn = self._GetReactingAtomsMN(1)
        return self.reacting_atoms_mn

    def Canonicalize(self, depth=5):
        """ Canonicalize the reaction, the atoms in self.rxn
        (RDKit chemical reaction object) will be reordered. Then the output
        smarts will be canonical, the canonicalization will be perfect as
        depth -> infinity.

        :param depth:
        :return:
        """
        can_agents = []
        for agent in self.agents:
            can_agents.append(agent)
        can_reactants = []
        for reactant in self.reactants:
            can_reactants.append(
                self._CanonicalizeMol(reactant, depth=depth))
            self._RemoveAtomMap(reactant)
        self.rxn.RemoveUnmappedReactantTemplates(thresholdUnmappedAtoms=1e-5)
        for reactant in can_reactants:
            self.rxn.AddReactantTemplate(reactant)
        can_products = []
        for product in self.products:
            can_products.append(
                self._CanonicalizeMol(product, depth=depth))
            self._RemoveAtomMap(product)
        self.rxn.RemoveUnmappedProductTemplates(thresholdUnmappedAtoms=1e-5)
        for product in can_products:
            self.rxn.AddProductTemplate(product)
        self.rxn.RemoveAgentTemplates()
        map(self.rxn.AddAgentTemplate, can_agents)
        return self.rxn

    def GetChiralityInfo(self):
        # get chirality_record of reactants and products.
        reactants_chirality_record = self._GetMolsChiralityRecord(
            self.reactants, 'reactant'
        )
        products_chirality_record = self._GetMolsChiralityRecord(
            self.products, 'product'
        )
        self.reactants_chirality_record = reactants_chirality_record
        self.products_chirality_record = products_chirality_record
        return reactants_chirality_record, products_chirality_record

    def RearangeMols(self):
        for reactant in self.reactants:
            if not self._IsReactMol(reactant, self.ReactingAtomsMN):
                self._RemoveAtomMap(reactant)
        for product in self.products:
            if not self._IsReactMol(product, self.ReactingAtomsMN):
                self._RemoveAtomMap(product)
        self.rxn.RemoveUnmappedReactantTemplates(thresholdUnmappedAtoms=1e-5)
        self.rxn.RemoveUnmappedProductTemplates(thresholdUnmappedAtoms=1e-5)

    def ReassignMappingNumber(self, depth=5):
        # Delete map number that exists only in reactants or products.
        reactantsAtomMapList = self._GetAtomMapList(self.reactants)
        productsAtomMapList = self._GetAtomMapList(self.products)
        for reactant in self.reactants:
            for atom in reactant.GetAtoms():
                AMN = atom.GetPropsAsDict().get('molAtomMapNumber')
                if AMN is not None and AMN not in productsAtomMapList:
                    atom.ClearProp('molAtomMapNumber')
        for product in self.products:
            for atom in product.GetAtoms():
                AMN = atom.GetPropsAsDict().get('molAtomMapNumber')
                if AMN is not None and AMN not in reactantsAtomMapList:
                    atom.ClearProp('molAtomMapNumber')
        # sort labeled atoms in reactant based on atomic environment
        atoms = []
        AEs = []
        for reactant in self.reactants:
            for atom in reactant.GetAtoms():
                map_number = atom.GetPropsAsDict().get('molAtomMapNumber')
                if map_number is not None:
                    AE = AtomEnvironment(reactant, atom, depth=depth,
                                         IsSanitized=False)
                    AEs.append(AE)
                    atoms.append(atom)
        sort_idx = np.argsort(AEs)
        atoms = np.asarray(atoms)[sort_idx]
        # relabel molAtomMapNumber from 1 to N based on atom order.
        atoms_map_change_dict = dict()
        label = 1
        for atom in atoms:
            map_number = atom.GetPropsAsDict()['molAtomMapNumber']
            assert (map_number not in atoms_map_change_dict)
            atoms_map_change_dict[map_number] = label
            atom.SetAtomMapNum(label)
            label += 1
        for product in self.products:
            for atom in product.GetAtoms():
                AMN = getAtomMapNumber(atom)
                if AMN in atoms_map_change_dict:
                    atom.SetAtomMapNum(atoms_map_change_dict[AMN])

    def Sanitized(self):
        # The effect of this line is not sure.
        rdChemReactions.SanitizeRxn(self.rxn)
        # Sanitize all molecules.
        for mol in self.reactants + self.products + self.agents:
            Chem.SanitizeMol(mol)

    def _AddChiralOnSmarts(self, smarts, chirality_record):
        for mp, record in chirality_record.items():
            seq = re.findall(':(%s)]' % '|'.join(record[1]), smarts)
            if self._IsChiralRecordChange(record[1], seq):
                smarts = smarts.replace(':%d]' % mp,
                                        ';%s:%d]' % (record[0].reverse, mp))
            else:
                smarts = smarts.replace(':%d]' % mp,
                                        ';%s:%d]' % (record[0].value, mp))
        return smarts

    def _CanonicalizeMol(self, mol, depth):
        """ Reorder the atoms sequence in mol, this is helpful to make the output
        SMILES or SMARTS string canonical.

        Parameters
        ----------
        mol: RDKit molecule object
        IsSanitized: Set False only for molecular fragments in reaction template.

        Returns
        -------
        Canonicalzed RDKit molecule object.
        """
        # Get sorted labeled atoms and unlabeled atoms
        labeled_atoms = []
        labeled_AEs = []
        unlabeled_atoms = []
        unlabeled_AEs = []
        for atom in mol.GetAtoms():
            map_number = atom.GetPropsAsDict().get('molAtomMapNumber')
            if map_number is None:
                AE = AtomEnvironment(mol, atom, depth=depth,
                                     IsSanitized=not self.is_template,
                                     order_by_labeling=self.is_template)
                unlabeled_AEs.append(AE)
                unlabeled_atoms.append(atom)
            else:
                AE = AtomEnvironment(mol, atom, depth=depth,
                                     IsSanitized=not self.is_template,
                                     order_by_labeling=self.is_template)
                labeled_AEs.append(AE)
                labeled_atoms.append(atom)
        labeled_atoms = list(np.asarray(labeled_atoms)[np.argsort(labeled_AEs)])
        unlabeled_atoms = list(
            np.asarray(unlabeled_atoms)[np.argsort(unlabeled_AEs)])
        N_atoms = len(labeled_atoms) + len(unlabeled_atoms)

        # create an atom-sorted molecule
        m = Chem.MolFromSmiles('.'.join(['O'] * N_atoms))
        mw = Chem.RWMol(m)
        idx_old2new = dict()
        idx_new2old = dict()
        idx = 0
        for atom in labeled_atoms + unlabeled_atoms:
            idx_old2new[atom.GetIdx()] = idx
            idx_new2old[idx] = atom.GetIdx()
            mw.ReplaceAtom(idx, atom)
            idx += 1
        bonds_idx = []
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonds_idx.append((idx_old2new[i], idx_old2new[j]))
        bonds_idx.sort()
        for i, j in bonds_idx:
            bond = mol.GetBondBetweenAtoms(idx_new2old[i], idx_new2old[j])
            mw.AddBond(i, j, bond.GetBondType())
        can_mol = mw.GetMol()
        if not self.is_template:
            Chem.SanitizeMol(can_mol)
        return can_mol

    def _DictAtomMap2AtomEnvi(self, mols, depth):
        AtomMapDict = dict()
        for mol in mols:
            for atom in mol.GetAtoms():
                AMN = atom.GetPropsAsDict().get('molAtomMapNumber')
                if AMN is not None:
                    if not self.is_template:
                        AtomMapDict[AMN] = AtomEnvironment(
                            mol, atom, depth=depth, order_by_labeling=True)
                    else:
                        AtomMapDict[AMN] = AtomEnvironment(
                            mol, atom, depth=depth, order_by_labeling=True,
                            IsSanitized=False)
        return AtomMapDict

    def _GetMolsChiralityRecord(self, mols, category):
        chirality_record = dict()
        for mol in mols:
            for atom in mol.GetAtoms():
                if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                    if atom.GetPropsAsDict().get('molAtomMapNumber') is not None:
                        acr = self._GetAtomChiralityRecord(atom, category)
                        if acr is None:
                            warnings.warn('Neighbors of chiral atom '
                                          f'{atom.GetSmarts()} in %s are not '
                                          'labeled, clear the ChiralTag.'
                                          % category)
                            warnings.warn(f'reaction smarts: {self.reaction_smarts}')
                            atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                        else:
                            chirality_record[acr[0]] = [acr[1], acr[2]]
        return chirality_record

    def _GetAtomChiralityRecord(self, atom, category):
        neighbors_mn = []
        for neighbor in atom.GetNeighbors():
            mn = neighbor.GetPropsAsDict().get('molAtomMapNumber')
            if mn is None:
                # Cannot record the atomic chirality since one of its neighbors
                # is not labeled.
                return None
            neighbors_mn.append(str(mn))
        if category == 'reactant':
            smarts_string = self.reaction_smarts.split('>')[0]
        else:
            smarts_string = self.reaction_smarts.split('>')[-1]
        atom_mp = atom.GetPropsAsDict().get('molAtomMapNumber')
        # there is a hydrogen
        if len(neighbors_mn) == 3:
            neighbors_mn.append(str(atom_mp))
        atom_smarts = re.search('\[[0-9A-Za-z@:;&]+:%d]' % atom_mp,
                                smarts_string)[0]
        if '@@' in atom_smarts:
            record1 = ChiralTag('@@')
        elif '@' in atom_smarts:
            record1 = ChiralTag('@')
        else:
            record1 = None
        record2 = re.findall(':(%s)]' % '|'.join(neighbors_mn), smarts_string)
        return atom_mp, record1, record2

    def _GetReactingAtomsMN(self, depth):
        """

        :param depth: atom neighbor search depth.
        :return: The map number of a list of atoms that participate the
            reaction.
        """
        ReactingAtoms = []
        reactantAtomMap = self._DictAtomMap2AtomEnvi(self.reactants,
                                                     depth=depth)
        productAtomMap = self._DictAtomMap2AtomEnvi(self.products,
                                                    depth=depth)
        for idx, AEr in reactantAtomMap.items():
            AEp = productAtomMap.get(idx)
            if AEp is None:
                continue
            atom_r = AEr.tree.all_nodes()[0].data
            atom_p = AEp.tree.all_nodes()[0].data
            if AEr != AEp or \
                    (not self.is_template and atom_r.GetTotalNumHs() != atom_p.GetTotalNumHs()) or \
                    (self.is_template and atom_r.GetNumExplicitHs() != atom_p.GetNumExplicitHs()) or \
                    atom_r.GetFormalCharge() != atom_p.GetFormalCharge():
                ReactingAtoms.append(idx)
        return ReactingAtoms

    @staticmethod
    def _GetAtomMapList(mols):
        AtomMapList = []
        for mol in mols:
            for atom in mol.GetAtoms():
                AMN = atom.GetPropsAsDict().get('molAtomMapNumber')
                if AMN is not None:
                    assert (AMN not in AtomMapList)
                    AtomMapList.append(AMN)
        return AtomMapList

    @staticmethod
    def _GetAtomFromMapNumber(mols, map_number, return_mol=False):
        map_number = int(map_number)
        for mol in mols:
            for atom in mol.GetAtoms():
                if atom.GetPropsAsDict().get('molAtomMapNumber') == map_number:
                    if return_mol:
                        return atom, mol
                    else:
                        return atom
        return None

    @staticmethod
    def _IsChiralRecordChange(record1, record2):
        N_permutation = 0
        for i in range(len(record1)):
            index = record2.index(record1[i])
            if i != index:
                record2[i], record2[index] = record2[index], record2[i]
                N_permutation += 1
        if N_permutation % 2 == 0:
            return False
        else:
            return True

    @staticmethod
    def _IsReactMol(mol, ReactingAtoms):
        for atom in mol.GetAtoms():
            if atom.GetPropsAsDict().get('molAtomMapNumber') in ReactingAtoms:
                return True
        else:
            return False

    @staticmethod
    def _Mols2CanonicalSmiles(mols):
        fragments = []
        for mol in mols:
            smiles = Chem.MolToSmiles(mol)
            if smiles == '[HH2-]':
                smiles = '[H-]'
            fragments.append(smiles)
        fragments.sort()
        return '.'.join(fragments)

    @staticmethod
    def _Mols2CanonicalSmarts(mols):
        fragments = []
        for mol in mols:
            fragments.append(Base._Mol2CanonicalSmarts(mol))
        fragments.sort()
        return '.'.join(fragments)

    @staticmethod
    def _Mol2CanonicalSmarts(mol):
        atoms_to_use = []
        symbols = []
        for atom in mol.GetAtoms():
            atoms_to_use.append(atom.GetIdx())
            symbol = Base.__GetAtomSmarts(atom)
            symbols.append(symbol)
        smarts = Chem.MolFragmentToSmiles(
            mol, atoms_to_use,
            atomSymbols=symbols,
            allHsExplicit=True,
            isomericSmiles=True,
            allBondsExplicit=False)
        return smarts

    def _Frag2CanonicalSmarts(self, mols):
        fragments = []
        for mol in mols:
            atoms_to_use = []
            symbols = []
            for atom in mol.GetAtoms():
                atoms_to_use.append(atom.GetIdx())
                symbols.append(atom.GetSmarts().replace('&', ';'))
            fragments.append('(' + Chem.MolFragmentToSmiles(
                mol, atoms_to_use,
                atomSymbols=symbols,
                allHsExplicit=True,
                isomericSmiles=True,
                allBondsExplicit=True) + ')')
        fragments.sort()
        return '.'.join(fragments)

    @staticmethod
    def __GetAtomSmarts(atom):
        symbol = atom.GetSmarts()
        rep = {
            '&+': '+',
            '&-': '-',
            '&H1': 'H',
            '&H': 'H',
            '&': ';'
        }
        rep = dict((re.escape(k), v)
                   for k, v in rep.items())
        pattern = re.compile("|".join(rep.keys()))
        symbol = pattern.sub(lambda m: rep[re.escape(m.group(0))], symbol)
        return symbol

    @staticmethod
    def _GetStrictAtomSmarts(atom):
        ''' Copy from Rdchiral.
        For an RDkit atom object, generate a SMARTS pattern that
        matches the atom as strictly as possible
        '''
        symbol = atom.GetSmarts()
        if atom.GetSymbol() == 'H':
            symbol = '[#1]'

        if '[' not in symbol:
            symbol = '[' + symbol + ']'

        if 'H' not in symbol:
            H_symbol = 'H{}'.format(atom.GetTotalNumHs())
            # Explicit number of hydrogens: include "H0" when no hydrogens present
            if ':' in symbol:  # stick H0 before label
                symbol = symbol.replace(':', ';{}:'.format(H_symbol))
            else:
                symbol = symbol.replace(']', ';{}]'.format(H_symbol))

        # Explicit degree
        if ':' in symbol:
            symbol = symbol.replace(':', ';D{}:'.format(atom.GetDegree()))
        else:
            symbol = symbol.replace(']', ';D{}]'.format(atom.GetDegree()))

        # Explicit formal charge
        if '+' not in symbol and '-' not in symbol:
            charge = atom.GetFormalCharge()
            charge_symbol = '+' if (charge >= 0) else '-'
            charge_symbol += '{}'.format(abs(charge))
            if ':' in symbol:
                symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
            else:
                symbol = symbol.replace(']', ';{}]'.format(charge_symbol))

        return symbol

    @staticmethod
    def _GetGeneralAtomSmarts(atom):
        '''Copy from Rdchiral.
        This function takes an RDKit atom and turns it into a wildcard
        using heuristic generalization rules. This function should be used
        when candidate atoms are used to extend the reaction core for higher
        generalizability
        '''
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
    def _RemoveAtomMap(mol):
        if mol.__class__ == list:
            for m in mol:
                for atom in m.GetAtoms():
                    atom.ClearProp('molAtomMapNumber')
        else:
            for atom in mol.GetAtoms():
                atom.ClearProp('molAtomMapNumber')


class ChiralTag:
    def __init__(self, value):
        assert (value in ['@', '@@'])
        self.value = value

    @property
    def reverse(self):
        if self.value == '@':
            return '@@'
        else:
            return '@'
