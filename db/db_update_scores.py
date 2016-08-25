import os

from db_declarative import Reference, Base, Fasta, Read, Score
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker

from config import *

engine = create_engine(db_name)

Base.metadata.bin = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()


def read_score_file(file, alignment_type):
    if os.path.exists(file):
        print(file)
        with open(file, 'r') as f:
            f.readline()
            current_fasta = session.query(Fasta.id, Fasta.reference_id).filter(Fasta.reference_id == ref) \
                .filter(Fasta.coverage == c).filter(Fasta.length == l).first()
            reads = session.query(Read.id, Read.fasta_id, Read.index).filter(Read.fasta_id == current_fasta.id)
            read_ids = [elem[0] for elem in reads.all()]
            read_index2id = {elem[2]: elem[0] for elem in reads.all()}
            scores = {(score.tail_id, score.head_id): score.id for score in
                      session.query(Score.id, Score.head_id, Score.tail_id).filter().all() if Score.tail_id in read_ids}
            for line in f:
                x, y, z, = line.split('\t')
                tail = read_index2id[int(x)]
                head = read_index2id[int(y)]
                if alignment_type == 'calign':
                    new_score = Score(tail_id=tail, head_id=head, calign_score=int(z))
                elif alignment_type == 'seqalign':
                    new_score = Score(tail_id=tail, head_id=head, seq_align_score=int(z))
                elif alignment_type == 'seq5align':
                    new_score = Score(tail_id=tail, head_id=head, seq5_align_score=int(z))
                elif alignment_type == 'seq20align':
                    new_score = Score(tail_id=tail, head_id=head, seq20_align_score=int(z))
                if scores.get((tail, head)):
                    new_score.id = scores.get((tail, head))
                session.merge(new_score)
        session.commit()


def read_score_files(file):
    read_score_file(file + '.calignscore', 'calign')
    read_score_file(file + '.seqscore', 'seqalign')
    read_score_file(file + '_banded5.seqscore', 'seq5align')
    read_score_file(file + '_banded20.seqscore', 'seq20align')


for ref in [1]:
    new_ref = Reference(id=ref, name=references[ref - 1])
    session.merge(new_ref)
    for c in [40]:
        for l in [400]:
            file_name = DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}'.format(ref, c, l)
            read_score_files(file_name)
session.commit()
