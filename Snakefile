from pathlib import Path

configfile: 'config.yaml'

outputs = Path(config['outdir'])

rule all:
