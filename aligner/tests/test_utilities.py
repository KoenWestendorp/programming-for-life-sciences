import aligner.utilities as utilities


def test_read_fasta():
    raw = """>sp|P0DTC3|AP3A_SARS2 ORF3a protein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=3a PE=1 SV=1
MDLFMRIFTIGTVTLKQGEIKDATPSDFVRATATIPIQASLPFGWLIVGVALLAVFQSAS
KIITLKKRWQLALSKGVHFVCNLLLLFVTVYSHLLLVAAGLEAPFLYLYALVYFLQSINF
VRIIMRLWLCWKCRSKNPLLYDANYFLCWHTNCYDYCIPYNSVTSSIVITSGDGTTSPIS
EHDYQIGGYTEKWESGVKDCVVLHSYFTSDYYQLYSTQLSTDTGVEHVTFFIYNKIVDEP
EEHVQIHTIDGSSGVVNPVMEPIYDEPTTTTSVPL
>sp|P0DTC5|VME1_SARS2 Membrane protein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=M PE=1 SV=1
MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLWLLWPV
TLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMWSFNPETNILL
NVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLPKEITVATSRTLSYYK
LGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALLVQ
"""
    entries = utilities.read_fasta(raw)
    assert len(entries) == 2, "The raw fasta string contains two entries."
    assert entries[0] == utilities.FastaEntry(
        "sp|P0DTC3|AP3A_SARS2 ORF3a protein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=3a PE=1 SV=1",
        "MDLFMRIFTIGTVTLKQGEIKDATPSDFVRATATIPIQASLPFGWLIVGVALLAVFQSASKIITLKKRWQLALSKGVHFVCNLLLLFVTVYSHLLLVAAGLEAPFLYLYALVYFLQSINFVRIIMRLWLCWKCRSKNPLLYDANYFLCWHTNCYDYCIPYNSVTSSIVITSGDGTTSPISEHDYQIGGYTEKWESGVKDCVVLHSYFTSDYYQLYSTQLSTDTGVEHVTFFIYNKIVDEPEEHVQIHTIDGSSGVVNPVMEPIYDEPTTTTSVPL"
    )


