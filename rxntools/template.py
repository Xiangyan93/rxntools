#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .base import *


class ReactionTemplate(Base):
    def __init__(self, template_smarts):
        super().__init__(reaction_smarts=template_smarts, is_template=True)
        self.GetChiralityInfo()
        # canonical
        self.ReassignMappingNumber(depth=5)
        self.Canonicalize(depth=5)

    # Check whether the output canonical smarts keep identical.
    def Check(self):
        assert (ReactionTemplate(self.Smarts).Smarts == self.Smarts)