import sys
from ..Fasta_ToolsBox import get_N_percent

fasta = sys.argv[1]
output = sys.argv[2]
sample = sys.argv[3]

get_N_percent(fasta, output, sample)
