#!/usr/bin/env python3

import os
import sys
from series_utils.parse_dmpcas import calc_dmpci_file_hash

print(calc_dmpci_file_hash(sys.argv[1]))
