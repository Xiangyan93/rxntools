#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .base import *


class ReactionTemplate(Base):
    def __init__(self, template_smarts):
        super().__init__(reaction_smarts=template_smarts, is_template=True)
        self.GetChiralityInfo()
        # canonical
        self.ReassignMappingNumber(depth=100)
        self.Canonicalize(depth=100)
