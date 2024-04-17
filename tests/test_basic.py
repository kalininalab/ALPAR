from sr_amr.binary_tables import check_contigs

def test_basic1():
    assert check_contigs("./test.fna")

