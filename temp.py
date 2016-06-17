# dictionary of all samples
populations = {}

# wild population, w/ accession codes and sample ID
populations['wildfrench'] = {}
populations['wildfrench']['SRR997319']='Avey36' # Aveyron
populations['wildfrench']['SRR997317']='Fos6'   # Fos-su-Mer
populations['wildfrench']['SRR997304']='Fos2'   # Fos-su-Mer
populations['wildfrench']['SRR997303']='Her65'  # Herauld
populations['wildfrench']['SRR997318']='Lan7'   # Lancon
populations['wildfrench']['SRR997316']='Lan8'   # Lancon
populations['wildfrench']['SRR997305']='Vau73'  # Vaucluse

# domestic population, w/ accession codes and sample ID
populations['domestic'] = {}
populations['domestic']['SRR997325']='BH23'      # Belgian hare
populations['domestic']['SRR997320']='FA801'     # Champagne d'argent
populations['domestic']['SRR997321']='AC100'     # Angora
populations['domestic']['SRR997327']='A93015'    # Angora
populations['domestic']['SRR997323']='FL920'     # French lop
populations['domestic']['SRR997326']='FG3'       # Flemish giant
populations['domestic']['SRR997324']='FG4'       # Flemish giant
populations['domestic']['SRR997322']='REX12'     # Rex


for population, samples in populations.iteritems():
    for sample, code in samples.iteritems():
        print sample, code