import os


def create_db(dazz_db_dir, dest_db, source_files):
    """

    :param dest_db:
    :param source_files: list of fasta files
    :return:
    """
    os.system("{}/fasta2DB {} -f {}".format(dazz_db_dir, dest_db, " ".join(source_files)))


def align(daligner_dir, dir_of_files, filename):
    os.chdir(dir_of_files)
    os.system("{}/daligner {}.db {}.db".format(daligner_dir, filename, filename))
    os.system("{}/LAsort *.las".format(daligner_dir))
    os.system("{}/LAmerge {}.las *.S.las".format(daligner_dir, filename))


def create_ovl(daligner_dir, filename):
    print("Create overlap file")
    os.system("{}/LAdump -cdo {}.db {}.las > {}.ovl ".format(daligner_dir, filename, filename, filename))
