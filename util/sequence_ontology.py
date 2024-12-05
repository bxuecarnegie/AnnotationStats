from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag


PROTEIN_CODING = "SO:0000010"
PROTEIN_CODING_GENE = "SO:0001217"


def read_so_obo(so_obo_path, parent_so):
    godag = get_godag(so_obo_path)
    gosubdag = GoSubDag(parent_so, godag)
    go_id, go_term = max(gosubdag.go2obj.items(), key=lambda t: t[1].depth)
    print(go_id, go_term.name)




if __name__ == '__main__':
    so_obo = "/Users/bxue/Documents/Carnegie/SourceData/GFFs/D.melanogaster/FB2023_06/so-simple.obo"
    read_so_obo(so_obo, PROTEIN_CODING_GENE)
