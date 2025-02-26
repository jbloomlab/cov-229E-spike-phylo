
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/sharari/miniforge3/envs/evo_jump/lib/python3.11/site-packages', '/home/sharari/.cache/snakemake/snakemake/source-cache/runtime-cache/tmp6pxhuzue/file/fh/fast/bloom_j/computational_notebooks/sharari/2025/HcoV_229E_phylo_analysis/Rules/../Scripts', '/fh/fast/bloom_j/computational_notebooks/sharari/2025/HcoV_229E_phylo_analysis/Rules/../Scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xa3\x1b\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c1Configure/Input_Data/229E_full_accession_list.txt\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x0eaccession_list\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x1eResults/spike/nucleotide.fasta\x94\x8c\x1aResults/spike/metadata.tsv\x94e}\x94(h\x0c}\x94(\x8c\x0ffasta_sequences\x94K\x00N\x86\x94\x8c\x08metadata\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh*h&h,h\'ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(Mz\rM`mK\x00]\x94(\x8c\x08JX503061\x94\x8c\x08KY674919\x94\x8c\x08KY073748\x94\x8c\x08KY073748\x94e]\x94(\x8c\x08DQ243972\x94\x8c\x08DQ243976\x94\x8c\x08DQ243977\x94\x8c\x08KM055556\x94\x8c\x08KY369909\x94e]\x94(\x8c\tNC_028752\x94\x8c\tNC_002645\x94\x8c\x08DQ243972\x94\x8c\x08DQ243976\x94\x8c\x08DQ243977\x94\x8c\x08KM055556\x94\x8c\x08KY369909\x94e\x8c\x02No\x94]\x94(\x8c\x06Israel\x94\x8c\x07Lebanon\x94\x8c\x06Jordan\x94\x8c\x05Syria\x94\x8c\x07Georgia\x94\x8c\x04Iraq\x94\x8c\x07Armenia\x94\x8c\x0cSaudi Arabia\x94\x8c\x06Kuwait\x94\x8c\x05Yemen\x94\x8c\x07Bahrain\x94\x8c\x05Qatar\x94\x8c\x04Iran\x94\x8c\x14United Arab Emirates\x94\x8c\x04Oman\x94\x8c\nUzbekistan\x94\x8c\x0bAfghanistan\x94\x8c\nKazakhstan\x94\x8c\x08Pakistan\x94\x8c\nTajikistan\x94\x8c\x08Maldives\x94\x8c\nKyrgyzstan\x94\x8c\x05India\x94\x8c\tSri Lanka\x94\x8c\x05Nepal\x94\x8c\x04Asia\x94\x8c\nBangladesh\x94\x8c\x07Myanmar\x94\x8c\x08Thailand\x94\x8c\x08Malaysia\x94\x8c\x08Mongolia\x94\x8c\x04Laos\x94\x8c\tSingapore\x94\x8c\x05China\x94\x8c\x08Cambodia\x94\x8c\x07Vietnam\x94\x8c\tHong Kong\x94\x8c\tIndonesia\x94\x8c\x06Brunei\x94\x8c\x06Taiwan\x94\x8c\x0bPhilippines\x94\x8c\x0bTimor-Leste\x94\x8c\x0bSouth Korea\x94\x8c\x05Japan\x94e]\x94(\x8c\x10French Polynesia\x94\x8c\x05Palau\x94\x8c\tAustralia\x94\x8c\x10Papua New Guinea\x94\x8c\x07Oceania\x94\x8c\nMicronesia\x94\x8c\x0fSolomon Islands\x94\x8c\x07Vanuatu\x94\x8c\x10Marshall Islands\x94\x8c\x0bNew Zealand\x94\x8c\x08Kiribati\x94\x8c\x04Fiji\x94e]\x94(\x8c\nCabo Verde\x94\x8c\x06Gambia\x94\x8c\rGuinea-Bissau\x94\x8c\x07Senegal\x94\x8c\x0cSierra Leone\x94\x8c\x06Guinea\x94\x8c\x07Liberia\x94\x8c\nMauritania\x94\x8c\x10C\xc3\x83\xc2\xb4te d\'Ivoire\x94\x8c\x07Morocco\x94\x8c\x04Mali\x94\x8c\x0cBurkina Faso\x94\x8c\x05Ghana\x94\x8c\x04Togo\x94\x8c\x05Benin\x94\x8c\x07Algeria\x94\x8c\x15Sao Tome and Principe\x94\x8c\x07Nigeria\x94\x8c\x07Tunisia\x94\x8c\x05Niger\x94\x8c\x11Equatorial Guinea\x94\x8c\x05Gabon\x94\x8c\x08Cameroon\x94\x8c\x15Republic of the Congo\x94\x8c\x05Libya\x94\x8c\x07Namibia\x94\x8c\x06Angola\x94\x8c\x04Chad\x94\x8c\x18Central African Republic\x94\x8c Democratic Republic of the Congo\x94\x8c\x06Africa\x94\x8c\x08Botswana\x94\x8c\x0cSouth Africa\x94\x8c\x06Zambia\x94\x8c\x07Lesotho\x94\x8c\x05Sudan\x94\x8c\x0bSouth Sudan\x94\x8c\x08Zimbabwe\x94\x8c\x07Burundi\x94\x8c\x06Rwanda\x94\x8c\x05Egypt\x94\x8c\x08Eswatini\x94\x8c\x06Uganda\x94\x8c\x06Malawi\x94\x8c\nMozambique\x94\x8c\x08Tanzania\x94\x8c\x05Kenya\x94\x8c\x08Ethiopia\x94\x8c\x08Djibouti\x94\x8c\x14Union of the Comoros\x94\x8c\nMadagascar\x94\x8c\x07Somalia\x94\x8c\nSeychelles\x94\x8c\tMauritius\x94e]\x94(\x8c\x07Iceland\x94\x8c\x08Portugal\x94\x8c\x07Ireland\x94\x8c\x05Spain\x94\x8c\x0eUnited Kingdom\x94\x8c\x07Andorra\x94\x8c\x06France\x94\x8c\x07Belgium\x94\x8c\x0bNetherlands\x94\x8c\nLuxembourg\x94\x8c\x06Monaco\x94\x8c\x0bSwitzerland\x94\x8c\x06Norway\x94\x8c\x07Denmark\x94\x8c\rLiechtenstein\x94\x8c\x07Germany\x94\x8c\x06Europe\x94\x8c\x05Italy\x94\x8c\x05Malta\x94\x8c\x08Slovenia\x94\x8c\x06Sweden\x94\x8c\x07Austria\x94\x8c\x07Croatia\x94\x8c\x0eCzech Republic\x94\x8c\x16Bosnia and Herzegovina\x94\x8c\x06Poland\x94\x8c\x08Slovakia\x94\x8c\nMontenegro\x94\x8c\x07Albania\x94\x8c\x07Hungary\x94\x8c\x06Serbia\x94\x8c\x06Kosovo\x94\x8c\x0fNorth Macedonia\x94\x8c\x06Greece\x94\x8c\tLithuania\x94\x8c\x07Romania\x94\x8c\x08Bulgaria\x94\x8c\x06Latvia\x94\x8c\x07Estonia\x94\x8c\x07Finland\x94\x8c\x07Belarus\x94\x8c\x07Moldova\x94\x8c\x07Ukraine\x94\x8c\x06Cyprus\x94\x8c\x06Turkey\x94\x8c\nAzerbaijan\x94\x8c\x06Russia\x94e]\x94(\x8c\tArgentina\x94\x8c\x07Uruguay\x94\x8c\x05Chile\x94\x8c\x08Paraguay\x94\x8c\x07Bolivia\x94\x8c\rSouth America\x94\x8c\x06Brazil\x94\x8c\x04Peru\x94\x8c\x07Ecuador\x94\x8c\x08Colombia\x94\x8c\x08Suriname\x94\x8c\x06Guyana\x94\x8c\tVenezuela\x94\x8c\x13Trinidad and Tobago\x94\x8c\x07Bonaire\x94\x8c\x07Curacao\x94\x8c\x05Aruba\x94e]\x94(\x8c\x06Panama\x94\x8c\nCosta Rica\x94\x8c\x07Grenada\x94\x8c\tNicaragua\x94\x8c Saint Vincent and the Grenadines\x94\x8c\x08Barbados\x94\x8c\x0bEl Salvador\x94\x8c\x0bSaint Lucia\x94\x8c\x08Honduras\x94\x8c\x08Dominica\x94\x8c\tGuatemala\x94\x8c\nGuadeloupe\x94\x8c\x06Belize\x94\x8c\x13Antigua and Barbuda\x94\x8c\x15Saint Kitts and Nevis\x94\x8c\x13Saint Barth\xc3\x83\xc2\xa9lemy\x94\x8c\x0cSint Maarten\x94\x8c\x0cSaint Martin\x94\x8c\x07Jamaica\x94\x8c\x12Dominican Republic\x94\x8c\x05Haiti\x94\x8c\x06Mexico\x94\x8c\x04Cuba\x94\x8c\x07Bahamas\x94\x8c\rNorth America\x94\x8c\x07Bermuda\x94\x8c\x03USA\x94\x8c\x06Canada\x94e]\x94(\x8c\x0cassembly_gap\x94\x8c\x08C_region\x94\x8c\x03CDS\x94\x8c\ncentromere\x94\x8c\x06D-loop\x94\x8c\tD_segment\x94\x8c\x04exon\x94\x8c\x03gap\x94\x8c\x04gene\x94\x8c\x04iDNA\x94\x8c\x06intron\x94\x8c\tJ_segment\x94\x8c\x0bmat_peptide\x94\x8c\x0cmisc_binding\x94\x8c\x0fmisc_difference\x94\x8c\x0cmisc_feature\x94\x8c\x0bmisc_recomb\x94\x8c\x08misc_RNA\x94\x8c\x0emisc_structure\x94\x8c\x0emobile_element\x94\x8c\rmodified_base\x94\x8c\x04mRNA\x94\x8c\x05ncRNA\x94\x8c\x08N_region\x94\x8c\x0cold_sequence\x94\x8c\x06operon\x94\x8c\x04oriT\x94\x8c\npolyA_site\x94\x8c\rprecursor_RNA\x94\x8c\x0fprim_transcript\x94\x8c\x0bprimer_bind\x94\x8c\npropeptide\x94\x8c\x0cprotein_bind\x94\x8c\nregulatory\x94\x8c\rrepeat_region\x94\x8c\nrep_origin\x94\x8c\x04rRNA\x94\x8c\x08S_region\x94\x8c\x0bsig_peptide\x94\x8c\x06source\x94\x8c\tstem_loop\x94\x8c\x03STS\x94\x8c\x08telomere\x94\x8c\x05tmRNA\x94\x8c\x0ftransit_peptide\x94\x8c\x04tRNA\x94\x8c\x06unsure\x94\x8c\x08V_region\x94\x8c\tV_segment\x94\x8c\tvariation\x94\x8c\x053\'UTR\x94\x8c\x055\'UTR\x94e}\x94\x8c\x05Human\x94]\x94(\x8c\x0chomo sapiens\x94\x8c\x06humans\x94\x8c\x05human\x94\x8c\x07patient\x94es]\x94(\x8c\x01S\x94\x8c\rspike protein\x94\x8c\x0eS glycoprotein\x94\x8c\x12Spike glycoprotein\x94\x8c\x12spike glycoprotein\x94\x8c\x05spike\x94\x8c\x0cglycoprotein\x94\x8c\x14surface glycoprotein\x94\x8c\x01G\x94\x8c\tG PROTEIN\x94\x8c\x17attachment glycoprotein\x94\x8c\tG protein\x94\x8c\x02HA\x94e\x8c\x05spike\x94e}\x94(h\x0c}\x94(\x8c\x1bgenome_size_threshold_lower\x94K\x00N\x86\x94\x8c\x1bgenome_size_threshold_upper\x94K\x01N\x86\x94\x8c\nmax_frac_N\x94K\x02N\x86\x94\x8c\x16accesstions_to_exclude\x94K\x03N\x86\x94\x8c\x19accessions_from_pre_study\x94K\x04N\x86\x94\x8c\x15accessions_to_include\x94K\x05N\x86\x94\x8c\x11remove_duplicates\x94K\x06N\x86\x94\x8c\x04asia\x94K\x07N\x86\x94\x8c\x07oceania\x94K\x08N\x86\x94\x8c\x06africa\x94K\tN\x86\x94\x8c\x06europe\x94K\nN\x86\x94\x8c\rsouth_america\x94K\x0bN\x86\x94\x8c\rnorth_america\x94K\x0cN\x86\x94\x8c\x10genbank_features\x94K\rN\x86\x94\x8c\x0chost_grouper\x94K\x0eN\x86\x94\x8c\ngene_names\x94K\x0fN\x86\x94\x8c\x0cdesired_gene\x94K\x10N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bjl\x01\x00\x00Mz\rjn\x01\x00\x00M`mjp\x01\x00\x00K\x00jr\x01\x00\x00h;jt\x01\x00\x00h@jv\x01\x00\x00hFjx\x01\x00\x00hNjz\x01\x00\x00hOj|\x01\x00\x00h|j~\x01\x00\x00h\x89j\x80\x01\x00\x00h\xc0j\x82\x01\x00\x00h\xf0j\x84\x01\x00\x00j\x02\x01\x00\x00j\x86\x01\x00\x00j\x1f\x01\x00\x00j\x88\x01\x00\x00jT\x01\x00\x00j\x8a\x01\x00\x00j[\x01\x00\x00j\x8c\x01\x00\x00ji\x01\x00\x00ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94ji\x01\x00\x00a}\x94(h\x0c}\x94\x8c\x04gene\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94b\x8c\x04gene\x94ji\x01\x00\x00ub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x15/loc/scratch/11600167\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bj\xb1\x01\x00\x00K\x01j\xb3\x01\x00\x00K\x01j\xb5\x01\x00\x00j\xae\x01\x00\x00ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8c6Results/Logs/spike/download_and_process_accessions.txt\x94a}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0eAccession_list\x94\x8c1Configure/Input_Data/229E_full_accession_list.txt\x94\x8c\x11Reference_genbank\x94\x8c4Configure/Input_Data/HcoV-229E_referance_sequence.gb\x94\x8c\x13Reference_accession\x94\x8c\tNC_002645\x94\x8c\x08Outgroup\x94\x8c,camel_Riyadh_Ry141_2015_NC_028752_2015-03-XX\x94\x8c\x1bGenome_size_threshold_lower\x94Mz\r\x8c\x1bGenome_size_threshold_upper\x94M`m\x8c\x0cORF_min_size\x94}\x94\x8c\x05spike\x94Mz\rs\x8c\nmax_frac_N\x94}\x94\x8c\x05genes\x94K\x00s\x8c\x11Remove_duplicates\x94hN\x8c\x15Accessions_to_exclude\x94]\x94(h<h=h>h?e\x8c\x15Accessions_to_include\x94]\x94(hGhHhIhJhKhLhMe\x8c\x19Accessions_from_pre_study\x94]\x94(hAhBhChDhEe\x8c\x0eGene_groupings\x94}\x94\x8c\x05spike\x94j[\x01\x00\x00s\x8c\x0eAuspice_config\x94\x8c(Configure/Input_Data/auspice_config.json\x94\x8c\rColor_schemes\x94\x8c&Configure/Input_Data/color_schemes.tsv\x94\x8c\x0fColor_orderings\x94\x8c(Configure/Input_Data/color_orderings.tsv\x94\x8c\tLat_longs\x94\x8c"Configure/Input_Data/lat_longs.tsv\x94\x8c\x13Auspice_tree_titles\x94}\x94\x8c\x05spike\x94\x8c\x1c\'hCoV 229E phylogentic tree\'\x94s\x8c\x19Auspice_tree_descriptions\x94}\x94\x8c\x05spike\x94\x8c)Configure/Input_Data/spike_description.md\x94s\x8c\x04Asia\x94]\x94(hPhQhRhShThUhVhWhXhYhZh[h\\h]h^h_h`hahbhchdhehfhghhhihjhkhlhmhnhohphqhrhshthuhvhwhxhyhzh{e\x8c\x07Oceania\x94]\x94(h}h~h\x7fh\x80h\x81h\x82h\x83h\x84h\x85h\x86h\x87h\x88e\x8c\x06Africa\x94]\x94(h\x8ah\x8bh\x8ch\x8dh\x8eh\x8fh\x90h\x91h\x92h\x93h\x94h\x95h\x96h\x97h\x98h\x99h\x9ah\x9bh\x9ch\x9dh\x9eh\x9fh\xa0h\xa1h\xa2h\xa3h\xa4h\xa5h\xa6h\xa7h\xa8h\xa9h\xaah\xabh\xach\xadh\xaeh\xafh\xb0h\xb1h\xb2h\xb3h\xb4h\xb5h\xb6h\xb7h\xb8h\xb9h\xbah\xbbh\xbch\xbdh\xbeh\xbfe\x8c\x06Europe\x94]\x94(h\xc1h\xc2h\xc3h\xc4h\xc5h\xc6h\xc7h\xc8h\xc9h\xcah\xcbh\xcch\xcdh\xceh\xcfh\xd0h\xd1h\xd2h\xd3h\xd4h\xd5h\xd6h\xd7h\xd8h\xd9h\xdah\xdbh\xdch\xddh\xdeh\xdfh\xe0h\xe1h\xe2h\xe3h\xe4h\xe5h\xe6h\xe7h\xe8h\xe9h\xeah\xebh\xech\xedh\xeeh\xefe\x8c\rSouth America\x94]\x94(h\xf1h\xf2h\xf3h\xf4h\xf5h\xf6h\xf7h\xf8h\xf9h\xfah\xfbh\xfch\xfdh\xfeh\xffj\x00\x01\x00\x00j\x01\x01\x00\x00e\x8c\rNorth America\x94]\x94(j\x03\x01\x00\x00j\x04\x01\x00\x00j\x05\x01\x00\x00j\x06\x01\x00\x00j\x07\x01\x00\x00j\x08\x01\x00\x00j\t\x01\x00\x00j\n\x01\x00\x00j\x0b\x01\x00\x00j\x0c\x01\x00\x00j\r\x01\x00\x00j\x0e\x01\x00\x00j\x0f\x01\x00\x00j\x10\x01\x00\x00j\x11\x01\x00\x00j\x12\x01\x00\x00j\x13\x01\x00\x00j\x14\x01\x00\x00j\x15\x01\x00\x00j\x16\x01\x00\x00j\x17\x01\x00\x00j\x18\x01\x00\x00j\x19\x01\x00\x00j\x1a\x01\x00\x00j\x1b\x01\x00\x00j\x1c\x01\x00\x00j\x1d\x01\x00\x00j\x1e\x01\x00\x00e\x8c\x10Genbank_features\x94]\x94(j \x01\x00\x00j!\x01\x00\x00j"\x01\x00\x00j#\x01\x00\x00j$\x01\x00\x00j%\x01\x00\x00j&\x01\x00\x00j\'\x01\x00\x00j(\x01\x00\x00j)\x01\x00\x00j*\x01\x00\x00j+\x01\x00\x00j,\x01\x00\x00j-\x01\x00\x00j.\x01\x00\x00j/\x01\x00\x00j0\x01\x00\x00j1\x01\x00\x00j2\x01\x00\x00j3\x01\x00\x00j4\x01\x00\x00j5\x01\x00\x00j6\x01\x00\x00j7\x01\x00\x00j8\x01\x00\x00j9\x01\x00\x00j:\x01\x00\x00j;\x01\x00\x00j<\x01\x00\x00j=\x01\x00\x00j>\x01\x00\x00j?\x01\x00\x00j@\x01\x00\x00jA\x01\x00\x00jB\x01\x00\x00jC\x01\x00\x00jD\x01\x00\x00jE\x01\x00\x00jF\x01\x00\x00jG\x01\x00\x00jH\x01\x00\x00jI\x01\x00\x00jJ\x01\x00\x00jK\x01\x00\x00jL\x01\x00\x00jM\x01\x00\x00jN\x01\x00\x00jO\x01\x00\x00jP\x01\x00\x00jQ\x01\x00\x00jR\x01\x00\x00jS\x01\x00\x00e\x8c\x0cHost_grouper\x94jT\x01\x00\x00u\x8c\x04rule\x94\x8c\x1fdownload_and_process_accessions\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c_/fh/fast/bloom_j/computational_notebooks/sharari/2025/HcoV_229E_phylo_analysis/Rules/../Scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/fh/fast/bloom_j/computational_notebooks/sharari/2025/HcoV_229E_phylo_analysis/Scripts/download_NCBI_sequences.py';
######## snakemake preamble end #########
# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr

# Imports
import re
import datetime
import pandas as pd
from Bio import Entrez
Entrez.email = "229Eexample@229Eexample.com"     # Always tell NCBI who you are

# Functions
def read_and_process_accession_list(
    input_file_name, 
    output_fasta_file_name, 
    output_metadata_file_name, 
    length_threshold, 
    remove_duplicates,
    desired_gene,
    output_log_file,
    ):
    """
    Function to read in list of accessions, download genbank files,
    parse genbank files, and extract sequences/metadata
    """

    def check_if_genbank_feature(feature):
        """
        Function to check if a feature is
        a genbank feature
        """

        if feature in snakemake.params.genbank_features:
            return True
        else:
            return False

    def country_extraction(location):
        """
        Function to extract country
        """
        # lower case location and extract before colon
        location = location.lower().split(":")[0]

        # Intialize results
        extracted_country = None
        region = None

        # Reformat specific countries
        if location == "viet nam":
            extracted_country = "Vietnam"
        elif location == "usa":
            extracted_country = "USA"
        elif location == "zaire":
            extracted_country = "Democratic Republic of the Congo"
        elif location == "democratic republic of the congo":
            extracted_country = "Democratic Republic of the Congo"
        elif location == "cote d'ivoire":
            extracted_country = "Cote d'Ivoire"
        elif location == "bosnia and herzegovina":
            extracted_country = "Bosnia and Herzegovina"
        elif location == "ussr":
            extracted_country = "Russia"
        else:
            # Capitalize first char in string
            formatted_list = []
            for substring in location.split(" "):
                formatted_list.append(substring.capitalize())
            
            extracted_country = " ".join(formatted_list)

        # Group by region
        if extracted_country in snakemake.params.asia:
            region = "Asia"
        elif extracted_country in snakemake.params.oceania:
            region = "Oceania"
        elif extracted_country in snakemake.params.africa:
            region = "Africa"
        elif extracted_country in snakemake.params.europe:
            region = "Europe"
        elif extracted_country in snakemake.params.south_america:
            region = "South America"
        elif extracted_country in snakemake.params.north_america:
            region = "North America"
        else:
            region = "?"
        
        # Return results
        return (region, extracted_country)
    
    
    def host_grouper(host):
        """
        Function to group host species
        """

        # Initialize dict and results
        host_dict = snakemake.params.host_grouper
        result = "UNKNOWN"

        # Iterate through host dict
        for grouped_host, host_list in host_dict.items():
            if host.lower() in host_list:
                result = grouped_host
                break

        # Return results
        return result
    

    def parse_CDS(raw_CDS_string, sequence):
        """
        Function to extract CDS
        """

        # Extract CDS coordinates
        temp = re.findall(r"\d+", raw_CDS_string)
        CDS_coordinates = list(map(int, temp))

        # Check to make sure there is an even number of coordinates
        assert len(CDS_coordinates) % 2 == 0, "Odd number of coordinates"

        # Intialize CDS string
        CDS_string = ""
        n = len(sequence)

        # Make reverse complement of sequence
        sequence_RC = sequence[::-1] 
        sequence_RC = sequence_RC.upper()
        sequence_RC = (
            sequence_RC
            .replace("A", "t")
            .replace("C", "g")
            .replace("T", "a")
            .replace("G", "c")
            .replace("N", "n")
        )

        # Process each CDS coordinates
        if "join" in raw_CDS_string and "complement" in raw_CDS_string:
            for i in range(len(CDS_coordinates)-1, -1, -2):
                CDS_string += sequence_RC[n - CDS_coordinates[i]:n - CDS_coordinates[i-1] + 1]
        elif "join" in raw_CDS_string:
            for i in range(0, len(CDS_coordinates)-1, 2):
                CDS_string += sequence[CDS_coordinates[i]-1:CDS_coordinates[i+1]]
        elif "complement" in raw_CDS_string:
            CDS_string += sequence_RC[n - CDS_coordinates[1]:n - CDS_coordinates[0] + 1]
        else:
            CDS_string += sequence[CDS_coordinates[0]-1:CDS_coordinates[1]]

        return CDS_string


     # Open output files
    output_fasta_file = open(output_fasta_file_name, "w")
    output_metadata_file = open(output_metadata_file_name, "w")
    log_file = open(output_log_file, "w")
    total_accessions_count = 0
    removed_accessions_count = 0
    seen_sequences = []

    # Write header for metadata file
    header = [
        "strain",
        "virus",
        "gene",
        "host",
        "accession",
        "date",
        "location",
        "region",
        "country",
        "database",
        "authors",
        "url",
        "title",
        "journal",
        "paper_url",
    ]
    header = "\t".join(header) + "\n"
    output_metadata_file.write(header)
    
    # Open file with list of accession numbers
    with open(input_file_name, "r") as input_file:
        for line in input_file:

            # Initialize results for metadata
            strain = "MISSING"
            virus = "?"
            gene = "?"
            host = "?"
            accession = "?"
            date = "?"
            location = "?"
            region = "?"
            country = "?"
            database = "?"
            authors = "?"
            url = "?"
            title = "?"
            journal = "?"
            paper_url = "?"
            features_flag = False
            sequence_flag = False
            nucleotide_sequence = ""
            curr_CDS = False
            CDS = ""
            reference_count = 0
            length = 0
            backup_date = ""
            previous_study = "no"

            # Extract current accession ID
            accession_from_list = line.split()[0]
            total_accessions_count += 1

            # Check if accession is to be excluded
            if accession_from_list[:-2] in snakemake.params.accesstions_to_exclude:
                log_file.write(f"{accession_from_list} excluded based on config file!\n")
                removed_accessions_count += 1
                continue

            # Retrieve genbank file for accession ID
            entrez_genbank = Entrez.efetch(
                db="nucleotide", 
                id=accession_from_list, 
                rettype="genbank", 
                retmode="text"
                )
            log_file.write("\n")
            log_file.write(f"Processing {accession_from_list}\n")
            # Parse genbank file line by line to retrieve all metadata
            for line in entrez_genbank:

                # Process line by removing spaces
                line = " ".join([ele for ele in line.split(" ") if ele != ""])
                split_line = line.split(" ")

                # Check if not in CDS feature 
                if curr_CDS and check_if_genbank_feature(split_line[0]):
                    curr_CDS = False

                # Extract current CDS if segment not found and 
                # CDS does not contain a complement sequence b/c rabies 
                # should be all from 5' direction
                if split_line[0] == "CDS" and gene != desired_gene:
                    CDS = line.split(" ")[1]
                    curr_CDS = True

                # Extract CDS products and determine which segment they come from
                if desired_gene != None and ("/product" in line or "/gene" in line) and curr_CDS and gene != desired_gene:
                    gene = line.replace("\n", "").replace("\"", "").split("=")[1].upper()
                    def check_list(string_list, word):
                        """
                        Function to check if a word has any
                        matching strings in a list
                        """
                        for substring in string_list:
                            if substring == word:
                                return True
                        return False

                    if check_list(snakemake.params.gene_names, gene):
                        log_file.write(f"{gene} is {desired_gene}!\n")
                        gene = desired_gene
                    else:
                        log_file.write(f"ERROR: {gene} product not known or not {desired_gene}\n")

                # Extract feature information from genbank file
                if features_flag == True:
                    if "/isolate" in line or "/strain" in line:
                        strain = line.replace("\n", "").replace("\"", "").split("=")[1]
                        # Remove special characters ()[]{}|#>< not desired in nextstrain
                        strain = strain.replace("@", "-")
                        for char in ["(",")","[","]","{","}","|","#",">","<",";"]:
                            strain = strain.replace(char, "")
                    if "/host" in line or "/lab_host" in line:
                        host = line.replace("\n", "").replace("\"", "").split("=")[1].split(";")[0]
                        grouped_host = host_grouper(host)
                        if grouped_host == "UNKNOWN":
                            log_file.write(f"{host} not unknown!\n")
                            host = "?"
                        else:
                            host = grouped_host
                    if "/geo_loc_name" in line:
                        # For now, just using a single location field
                        location = line.replace("\n", "").replace("\"", "").split("=")[1]
                        region_and_country = country_extraction(location)
                        region = region_and_country[0]
                        country = region_and_country[1]
                    if "/collection_date" in line:
                        unformatted_date = line.replace("\n", "").replace("\"", "").split("=")[1]
                        if len(unformatted_date.split("-")) == 1:
                            if "/" in unformatted_date:
                                date = datetime.datetime.strptime(unformatted_date.split("/")[0], "%Y").strftime("%Y") + "-XX-XX"
                            else:
                                date = datetime.datetime.strptime(unformatted_date, "%Y").strftime("%Y") + "-XX-XX"
                        elif len(unformatted_date.split("-")) == 2:
                            for fmt in ["%b-%Y", "%Y-%m"]:
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                except:
                                    log_file.write(f"ERROR! {unformatted_date} not recongized as {fmt}\n")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m") + "-XX"
                                    break
                        elif len(unformatted_date.split("-")) == 3:
                            for fmt in ["%d-%b-%Y", "%Y-%m-%d"]:
                                try:
                                    datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                except:
                                    log_file.write(f"ERROR! {unformatted_date} not recongized as {fmt}\n")
                                    continue
                                else:
                                    date = datetime.datetime.strptime(unformatted_date, fmt).strftime("%Y-%m-%d")
                                    break
                        else:
                            log_file.write("Datetime not between 1 and 3\n")

                # Extract nucleotide sequence
                if sequence_flag == True and split_line[0] != "//": 
                    nucleotide_sequence += "".join(split_line[1:]).replace("\n", "")
                
                if split_line[0] == "LOCUS":
                    # Create a backup date based on top line of genbank
                    unformatted_backup_date = split_line[-1].replace("\n", "")
                    if len(unformatted_backup_date.split("-")) == 1:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%Y").strftime("%Y") + "-XX-XX"
                    elif len(unformatted_backup_date.split("-")) == 2:
                        backup_date = datetime.datetime.strptime(unformatted_backup_date, "%b-%Y").strftime("%Y-%m") + "-XX"
                    elif len(unformatted_backup_date.split("-")) == 3:
                        for fmt in ["%d-%b-%Y", "%Y-%m-%d"]:
                            try:
                                datetime.datetime.strptime(unformatted_backup_date, fmt).strftime("%Y-%m-%d")
                            except:
                                log_file.write(f"ERROR! {unformatted_backup_date} not recongized as {fmt}\n")
                                continue
                            else:
                                backup_date = datetime.datetime.strptime(unformatted_backup_date, fmt).strftime("%Y-%m-%d")
                                break
                    else:
                        log_file.write("Datetime not between 1 and 3\n")
                    # Get sequence length
                    length = int(split_line[2])
                    continue
                if split_line[0] == "DEFINITION":
                    if any(word in split_line[1] for word in snakemake.params.host_grouper):
                        host = 'Human'
                    else:
                        continue
                if split_line[0] == "ACCESSION":
                    accession = split_line[1].replace("\n", "")
                    genbank_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
                    url = genbank_base_url + accession
                    continue
                if split_line[0] == "DBLINK":
                    database = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "REFERENCE":
                    reference_count += 1
                    continue
                if split_line[0] == "AUTHORS" and reference_count == 1:
                    authors = split_line[1].split(",")[0] + " et al"
                    authors = authors.replace("\n", "")
                    continue
                if split_line[0] == "TITLE" and reference_count == 1:
                    title = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "JOURNAL" and reference_count == 1:
                    journal = " ".join(split_line[1:]).replace("\n", "")
                    continue
                if split_line[0] == "PUBMED" and reference_count == 1:
                    pubmed_base_url = "https://pubmed.ncbi.nlm.nih.gov/"
                    paper_url = pubmed_base_url + split_line[1].replace("\n", "")
                    continue
                if split_line[0] == "FEATURES":
                    features_flag = True # Set features flag to then extract feature information
                    continue 
                if split_line[0] == "ORIGIN":
                    sequence_flag = True # Set sequence flag to extract nucleotide sequence
                    continue
                

            # Check length of sequence and do not add if below threshold
            if length < length_threshold[0] or length > length_threshold[1]:
                log_file.write(f"{accession_from_list} excluded because {length} not within length limits!\n")
                removed_accessions_count += 1
                continue

            # Check if only specific genes should be kept
            if desired_gene != None and gene != desired_gene:
                log_file.write(f"{accession_from_list} excluded because {gene} is not the desired gene ({desired_gene})!\n")
                removed_accessions_count += 1
                continue

            # Check if host is not human
            if host == "?" and accession_from_list[:-2] not in snakemake.params.accessions_to_include:
                log_file.write(f"{accession_from_list} excluded because non human host\n")
                removed_accessions_count += 1
                continue

            # Check if host is not human
            if accession_from_list[:-2] in snakemake.params.accessions_from_pre_study:
                previous_study = 'yes'
                continue

            # Check if CDS was found else extract CDS from full sequence
            if gene == desired_gene and CDS == "":
                log_file.write(f"{accession_from_list} excluded because no CDS was found!\n")
                removed_accessions_count += 1
                continue
            else:
                # Extract CDS sequence
                if desired_gene != None:
                    log_file.write(f"CDS found: {CDS}\n")
                    nucleotide_sequence = parse_CDS(CDS, nucleotide_sequence)

                if remove_duplicates == "Yes":
                    # Check if nulceotide sequence is a duplicate
                    if nucleotide_sequence in seen_sequences and accession_from_list[:-2] not in snakemake.params.accessions_to_include:
                        log_file.write(f"{accession_from_list} excluded because duplicate sequence!\n")
                        removed_accessions_count += 1
                        continue
                    else:
                        seen_sequences.append(nucleotide_sequence)

            # If not collection date is found, add date from top of genbank file
            if date == "?":
                date = backup_date

            # Check if sequence is below ambiguous base threshold
            if nucleotide_sequence.upper().count("N")/len(nucleotide_sequence) > snakemake.params.max_frac_N:
                log_file.write(f"{accession_from_list} excluded because of high N count!\n")
                removed_accessions_count += 1
                continue
            
            # Join strain/isolate name with accession and date to make sure it is unique
            if strain != "MISSING":
                strain = strain + "_" + accession + "_" + date
            else:
                strain = accession + "_" + date

            # Replace slashes, periods, and spaces in name with underscores
            strain = strain.replace("/", "_")
            strain = strain.replace(". ", "-")
            strain = strain.replace(" ", "_")
            strain = strain.replace(".", "-")

            # Make sure every sequence has a fasta strain name
            assert strain != "MISSING", "Virus strain name is missing"

            # Create new metadata line
            new_metadata_line = "\t".join([
                strain,
                virus,
                gene,
                host,
                accession,
                date,
                location,
                region,
                country,
                database,
                authors,
                url,
                title,
                journal,
                paper_url,
                previous_study,
            ])
            new_metadata_line += "\n"

            # Write new metadata line
            output_metadata_file.write(new_metadata_line)

            # Write current fasta sequence to output file
            output_fasta_file.write(f">{strain}\n")
            output_fasta_file.write(f"{nucleotide_sequence}\n")

    log_file.write("\n")
    log_file.write(f"A total of {total_accessions_count} were processed and ")
    log_file.write(f"{total_accessions_count-removed_accessions_count} were retained!\n")   
    # Close files
    input_file.close()
    output_fasta_file.close()
    output_metadata_file.close()
    log_file.close()


def main():
    """
    Main method
    """

    # Input files
    list_of_accessions = str(snakemake.input.accession_list) 
    # Params
    length_threshold = (
        int(str(snakemake.params.genome_size_threshold_lower)), 
        int(str(snakemake.params.genome_size_threshold_upper))
    )
    # Empty String to None Conversion
    None_conversion = lambda i : i or None
    desired_gene = None_conversion(str(snakemake.params.desired_gene))
    remove_duplicates = str(snakemake.params.remove_duplicates)

    # Output files
    fasta_output = str(snakemake.output.fasta_sequences) 
    metadata_output = str(snakemake.output.metadata)
    output_log_file = str(snakemake.log)

    read_and_process_accession_list(
        list_of_accessions, 
        fasta_output, 
        metadata_output, 
        length_threshold, 
        remove_duplicates,
        desired_gene,
        output_log_file,
        )


if __name__ == "__main__":
    main()