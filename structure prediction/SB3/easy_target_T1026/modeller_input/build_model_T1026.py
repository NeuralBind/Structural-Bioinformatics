#!/usr/bin/env python

from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='alignT1.ali',
              knowns='6F2S', sequence='UKNP',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
