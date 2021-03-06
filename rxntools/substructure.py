#!/usr/bin/env python
# -*- coding: utf-8 -*-
from treelib import Tree
from .smiles import *


def getAtomMapNumber(atom):
    return atom.GetPropsAsDict().get('molAtomMapNumber') or -1


class MolecularTree:
    def __init__(self, order_by_labeling=False, IsSanitized=True):
        self.order_by_labeling = order_by_labeling
        self.IsSanitized = IsSanitized

    def get_tag(self, atom, order):
        if self.IsSanitized:
            tag = [atom.GetAtomicNum(),
                   atom.GetFormalCharge(),
                   atom.GetTotalNumHs(),
                   atom.GetIsAromatic(),
                   order]
        else:
            tag = [atom.GetAtomicNum(),
                   get_Charge_from_atom_smarts(atom.GetSmarts()),
                   get_Hcount_from_atom_smarts(atom.GetSmarts()),
                   atom.GetIsAromatic(),
                   order]
        if self.order_by_labeling:
            tag.append(getAtomMapNumber(atom))
        return tag

    def tree_grow(self, mol, tree, depth=5):
        for _ in range(depth):
            for node in tree.all_nodes():
                if not node.is_leaf():
                    continue
                for atom in node.data.GetNeighbors():
                    tree_id = tree._identifier
                    if atom.GetIdx() != node.predecessor(tree_id=tree_id):
                        order = mol.GetBondBetweenAtoms(
                            atom.GetIdx(),
                            node.data.GetIdx()
                        ).GetBondTypeAsDouble()
                        identifier = atom.GetIdx()
                        while tree.get_node(identifier) is not None:
                            identifier += len(mol.GetAtoms())
                        tree.create_node(
                            tag=self.get_tag(atom, order),
                            identifier=identifier,
                            data=atom,
                            parent=node.identifier
                        )
        return tree

    def __eq__(self, other):
        if self.get_rank_list() == other.get_rank_list():
            return True
        else:
            return False

    def __lt__(self, other):
        if self.get_rank_list() < other.get_rank_list():
            return True
        else:
            return False

    def __gt__(self, other):
        if self.get_rank_list() > other.get_rank_list():
            return True
        else:
            return False

    def get_rank_list(self):
        rank_list = []
        expand_tree = self.tree.expand_tree(mode=Tree.WIDTH, reverse=True)
        for identifier in expand_tree:
            rank_list += self.tree.get_node(identifier).tag
        return rank_list


class FunctionalGroup(MolecularTree):
    """Functional Group.

    atom0 -> atom1 define a directed bond in the molecule. Then the bond is
    removed and the functional group is defined as a multitree. atom1 is the
    root node of the tree.

    Parameters
    ----------
    mol : molecule object in RDKit

    atom0, atom1 : atom object in RDKit

    depth: the depth of the multitree.

    Attributes
    ----------
    tree : multitree represent the functional group
        each node has 3 important attributes: tag: [atomic number, bond order
        with its parent], identifier: atom index defined in RDKit molecule
        object, data: RDKit atom object.

    """

    def __init__(self, mol, atom0, atom1, depth=5, order_by_labeling=False,
                 IsSanitized=True):
        super().__init__(order_by_labeling, IsSanitized)
        self.mol = mol
        tree = Tree()
        tree.create_node(
            tag=self.get_tag(atom0, 0.),
            identifier=atom0.GetIdx(),
            data=atom0
        )
        order = mol.GetBondBetweenAtoms(
            atom0.GetIdx(),
            atom1.GetIdx()
        ).GetBondTypeAsDouble()
        tree.create_node(
            tag=self.get_tag(atom1, order),
            identifier=atom1.GetIdx(),
            data=atom1,
            parent=atom0.GetIdx()
        )
        self.tree = self.tree_grow(mol, tree, depth=depth)

    def get_bonds_list(self):
        bonds_list = []
        expand_tree = self.tree.expand_tree(mode=Tree.WIDTH, reverse=True)
        for identifier in expand_tree:
            i = identifier
            j = self.tree.get_node(identifier).predecessor(
                tree_id=self.tree._identifier)
            if j is None:
                continue
            ij = (min(i, j), max(i, j))
            bonds_list.append(ij)
        return bonds_list


class AtomEnvironment(MolecularTree):
    def __init__(self, mol, atom, depth=1, order_by_labeling=False,
                 IsSanitized=True):
        super().__init__(order_by_labeling, IsSanitized)
        self.mol = mol
        tree = Tree()
        tree.create_node(
            tag=self.get_tag(atom, 0.),
            identifier=atom.GetIdx(),
            data=atom
        )
        self.tree = self.tree_grow(mol, tree, depth=depth)

    def get_nth_neighbors(self, n=1):
        assert (n.__class__ == int and n >= 1)
        neighbors = []
        for node in self.tree.all_nodes():
            depth = self.tree.depth(node)
            if depth == n:
                neighbors.append(node.data)
        return neighbors


class UnmappedFragment(MolecularTree):
    def __init__(self, mol, atom0, atom1, order_by_labeling=False,
                 IsSanitized=True):
        super().__init__(order_by_labeling, IsSanitized)
        self.mol = mol
        tree = Tree()
        tree.create_node(
            tag=self.get_tag(atom0, 0.),
            identifier=atom0.GetIdx(),
            data=atom0
        )
        order = mol.GetBondBetweenAtoms(
            atom0.GetIdx(),
            atom1.GetIdx()
        ).GetBondTypeAsDouble()
        tree.create_node(
            tag=self.get_tag(atom1, order),
            identifier=atom1.GetIdx(),
            data=atom1,
            parent=atom0.GetIdx()
        )
        self.tree = self.tree_grow(mol, tree)

    def tree_grow(self, mol, tree):
        depth = 10
        for _ in range(depth):
            for node in tree.all_nodes():
                if not node.is_leaf():
                    continue
                for atom in node.data.GetNeighbors():
                    if getAtomMapNumber(atom) != -1:
                        continue
                    tree_id = tree._identifier
                    if atom.GetIdx() != node.predecessor(tree_id=tree_id):
                        order = mol.GetBondBetweenAtoms(
                            atom.GetIdx(),
                            node.data.GetIdx()
                        ).GetBondTypeAsDouble()
                        identifier = atom.GetIdx()
                        while tree.get_node(identifier) is not None:
                            identifier += len(mol.GetAtoms())
                        tree.create_node(
                            tag=self.get_tag(atom, order),
                            identifier=identifier,
                            data=atom,
                            parent=node.identifier
                        )
        if tree.depth() > depth:
            raise RuntimeError('It seems you are dealing with a reaction with '
                               'a huge unmapped fragment that connected to the '
                               'reaction center.')
        return tree

    def get_atoms_idx(self):
        atoms_idx = []
        expand_tree = self.tree.expand_tree(mode=Tree.WIDTH, reverse=True)
        for identifier in expand_tree:
            if not self.tree.get_node(identifier).is_root():
                atoms_idx.append(identifier)
        return atoms_idx
