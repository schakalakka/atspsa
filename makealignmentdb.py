import atspsa

import os


def create_db(dest, source_files):
    """

    :param dest:
    :param source_files: list of fasta files
    :return:
    """
    os.system("{}/fasta2DB {} -f {}".format(atspsa.DAZZ_DB, dest, " ".join(source_files)))


def align(db):
    filename = db.split(".db")[0]
    os.system("{}/daligner {} {}".format(atspsa.DALIGNER, db, db))
    os.system("{}/LAsort *.las".format(atspsa.DALIGNER))
    os.system("{}/LAmerge {}.las *.S.las".format(atspsa.DALIGNER, filename))


def create_ovl(db, las):
    filename = db.split(".db")[0]
    os.system("{}/LAdump -cdo {} {} > {} ".format(atspsa.DALIGNER, db, las, filename))
    return 0
