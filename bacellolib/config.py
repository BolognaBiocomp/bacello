BLASTDB="/seqdb/uniprot_sprot.fasta"

BLASTALPH = "ARNDCQEGHILKMFPSTWYV"
HSSPALPH = "VLIMFWYGAPSTCHRKQEND"

locmap = {"Cytoplasm": ("Cytoplasm", "GO:0005737"),
          "Nucleus": ("Nucleus", "GO:0005634"),
          "Secretory": ("Extracellular", "GO:0005615"),
          "Mitochondrion": ("Mitochondrion", "GO:0005739"),
          "Chloroplast": ("Chloroplast", "GO:0009507")}

GOINFO = {"GO:0005737": {"uniprot": {"location": {"value": "Cytoplasm"}},
                       "GO": {"type": "GO", "id": "GO:0005737",
                               "properties": {"term": "C:cytoplasm"},
                               "evidences": [{"code": "ECO:0000256|SAM:BaCelLo"}]}},
          "GO:0005634": {"uniprot": {"location": {"value": "Nucleus"}},
                       "GO": {"type": "GO", "id": "GO:0005634",
                              "properties": {"term": "C:nucleus"},
                              "evidences": [{"code": "ECO:0000256|SAM:BaCelLo"}]}},
          "GO:0005615": {"uniprot": {"location": {"value": "Secreted"}},
                       "GO": {"type": "GO", "id": "GO:0005615",
                              "properties": {"term": "C:extracellular space"},
                              "evidences": [{"code": "ECO:0000256|SAM:BaCelLo"}]}},
          "GO:0005739": {"uniprot": {"location": {"value": "Mitochondrion"}},
                       "GO": {"type": "GO", "id": "GO:0005739",
                              "properties": {"term": "C:mitochondrion"},
                              "evidences": [{"code": "ECO:0000256|SAM:BaCelLo"}]}},
          "GO:0009507": {"uniprot": {"location": {"value": "Chloroplast"}},
                       "GO": {"type": "GO", "id": "GO:0009507",
                              "properties": {"term": "C:chloroplast"},
                              "evidences": [{"code": "ECO:0000256|SAM:BaCelLo"}]}}}
