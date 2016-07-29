import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from db_declarative import Reference, Base, Fasta, Read, Score
from sqlalchemy import create_engine
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
            filtered_reads = session.query(Read, Fasta).filter(Fasta.reference_id == ref) \
                .filter(Fasta.coverage == c).filter(Fasta.length == l)
            for line in f.readlines():
                x, y, z, = line.split('\t')
                old_tail = filtered_reads.filter(Read.index == int(x)).first().Read.id
                old_head = filtered_reads.filter(Read.index == int(y)).first().Read.id
                old_score = session.query(Score).filter_by(tail_id=old_tail, head_id=old_head).first()
                if alignment_type == 'calign':
                    new_score = Score(tail_id=old_tail, head_id=old_head, calign_score=int(z))
                elif alignment_type == 'seqalign':
                    new_score = Score(tail_id=old_tail, head_id=old_head, seq_align_score=int(z))
                elif alignment_type == 'seq5align':
                    new_score = Score(tail_id=old_tail, head_id=old_head, seq5_align_score=int(z))
                elif alignment_type == 'seq20align':
                    new_score = Score(tail_id=old_tail, head_id=old_head, seq20_align_score=int(z))
                if old_score:
                    new_score.id = old_score.id
                session.merge(new_score)
        session.commit()


def read_score_files(file):
    read_score_file(file + '.calignscore', 'calign')
    read_score_file(file + '.seqscore', 'seqalign')
    read_score_file(file + '_banded5.seqscore', 'seq5align')
    read_score_file(file + '_banded20.seqscore', 'seq20align')


def read_time_file(file):
    if os.path.exists(file):
        with open(file, 'r') as f:
            return float(f.read())
    else:
        return None


def read_time_files(file):
    calign_time = read_time_file(file + '.cact')
    seq_align_time = read_time_file(file + '.act')
    seq5_align_time = read_time_file(file + '_banded5.act')
    seq20_align_time = read_time_file(file + '_banded20.act')
    return calign_time, seq_align_time, seq5_align_time, seq20_align_time


for ref in [1, 2, 3]:
    new_ref = Reference(id=ref, name=references[ref - 1])
    session.merge(new_ref)
    for c in coverages:
        for l in average_length_list:
            file_name = DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}'.format(ref, c, l)
            calign_time, seq_align_time, seq5_align_time, seq20_align_time = read_time_files(file_name)
            old_fasta = session.query(Fasta).filter(Fasta.reference_id == ref).filter(Fasta.coverage == c).filter(
                Fasta.length == l).first()
            new_fasta = Fasta(coverage=c, length=l, reference=new_ref, filepath=file_name, calign_time=calign_time,
                              seq_align_time=seq_align_time, seq5_align_time=seq5_align_time,
                              seq20_align_time=seq20_align_time)
            if old_fasta:
                new_fasta.id = old_fasta.id
            session.merge(new_fasta)
            session.commit()

            current_fasta = session.query(Fasta).filter(Fasta.reference_id == ref).filter(Fasta.coverage == c).filter(
                Fasta.length == l).first()
            reads = [(seqRecord.name, str(seqRecord.seq)) for seqRecord in
                     SeqIO.parse(file_name + '.fasta', "fasta", generic_dna)]
            for i, read in enumerate(reads):
                read_begin, read_end = read[0].split('/')[-1].split('_')
                read_length = int(read_end) - int(read_begin)
                old_read = session.query(Read).filter(Read.fasta_id == new_fasta.id).filter(Read.index == i).first()
                new_read = Read(sequence=read[1], length=read_length, begin=int(read_begin), end=int(read_end),
                                index=i, fasta_id=current_fasta.id)
                if old_read:
                    new_read.id = old_read.id
                session.merge(new_read)

                # read_score_files(file_name)
session.commit()
