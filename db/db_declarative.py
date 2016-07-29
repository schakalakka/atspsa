from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from config import *

Base = declarative_base()


class Reference(Base):
    __tablename__ = 'reference'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    fastas = relationship("Fasta", backref='reference')


class Fasta(Base):
    __tablename__ = 'fasta'
    id = Column(Integer, primary_key=True)
    reference_id = Column(Integer, ForeignKey('reference.id'))
    coverage = Column(Integer)
    length = Column(Integer)
    filepath = Column(String)
    calign_time = Column(Float)
    seq_align_time = Column(Float)
    seq20_align_time = Column(Float)
    seq5_align_time = Column(Float)
    reads = relationship('Read', backref='fasta')


class Read(Base):
    __tablename__ = 'read'
    id = Column(Integer, primary_key=True)
    fasta_id = Column(Integer, ForeignKey('fasta.id'))
    index = Column(Integer)
    sequence = Column(String)
    length = Column(Integer)
    begin = Column(Integer)
    end = Column(Integer)


class Score(Base):
    __tablename__ = 'score'
    id = Column(Integer, primary_key=True)
    head_id = Column(Integer, ForeignKey('read.id'))
    tail_id = Column(Integer, ForeignKey('read.id'))
    calign_score = Column(Integer)
    seq_align_score = Column(Integer)
    seq5_align_score = Column(Integer)
    seq20_align_score = Column(Integer)


engine = create_engine(db_name)

Base.metadata.create_all(engine)
