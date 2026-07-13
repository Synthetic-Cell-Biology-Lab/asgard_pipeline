lts = """EGFIHDBD_02570

MNACDEKI_00316
CCCLDPND_00768

EFIMHAMI_01168
PCFOPMMF_00766
LCJGLNGN_01126
IELMKNIH_01375
BOOEDNDH_02189
EDAGINKG_03361
NHHKBCLI_00896
KNFGCEEK_01559
JILBHOBO_01145
BNMOIDIA_00210
OGFGKDNL_02744
GAGJKOJD_02280
GDFAMCLG_02267

HDGNKJKO_00173
CIDNLPIE_01958
JEBHEBKO_00261

AEGAMCKF_02726
OJHPMOMD_01135
GJEANDHH_01317
ABKEDPFO_00285
CBCCFIPL_02104
OIOBNGCI_01127
HLNFHMNH_03600
AHBDPKAE_01822
KKNAGKOM_03519
NAJGBKGL_02756
MNGIOEKN_00485
KPLPGFEG_03910
EDBAHGEC_00528
KKMFLBAM_01450
EHFMEEDL_00686
LDGAGGOM_02587
MCONNPJD_01143
ABGDKFDE_01955
FEDPECBA_02910
NFECLIOE_01712

JMCHHPCH_01551
DIOOMLCJ_01489
HDAPCNDD_02123

NDAMOAGK_00476

ELCOMMCJ_00875

DPJHMCOC_00756
DEAKNLCD_01662
JCCFDLHG_01361
MGOGIFIP_02592
IOJDHAAC_00275
NDIMGIAA_00191
JBPHBBHP_00636
JNHFDLAO_01254
JLNMILLL_00759
CLCCMGGL_00105
FBKLGEOG_01354
GNNHGECN_01985"""


from Bio import SeqIO

wanted = {lt.strip() for lt in lts.splitlines() if lt.strip()}

found = set()

with open(
    "/home/anirudh/asgard_pipeline/database/protein_sets/ftsz/ftsz1_msa/selected_sequences.faa",
    "w",
) as out:
    for record in SeqIO.parse(
        "/home/anirudh/asgard_pipeline/database/collated/Version1/filtered/85comp10con/fasta/v1_cp85_con10.fasta",
        "fasta",
    ):
        if record.id in wanted:
            SeqIO.write(record, out, "fasta")
            found.add(record.id)

print(f"Found {len(found)}/{len(wanted)} locus tags")

missing = wanted - found
if missing:
    print("\nMissing locus tags:")
    for lt in sorted(missing):
        print(lt)
