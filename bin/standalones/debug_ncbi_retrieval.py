from Bio import Entrez

Entrez.email = "anirudhbantwalbaliga@gmail.com"
Entrez.api_key = "eed94abceca160fc601883bc7e508f604608"  # https://www.ncbi.nlm.nih.gov/account/

test_id = "MBM04725.1"

handle = Entrez.efetch(
    db="protein",
    id=test_id,
    retmode="xml"
)

data = Entrez.read(handle)
handle.close()

print(data)