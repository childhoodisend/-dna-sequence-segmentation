import gffutils

def parse():
    db = gffutils.create_db('source/sequence.gff3', 'ecoli.db', merge_strategy='merge', force=True)
    return db