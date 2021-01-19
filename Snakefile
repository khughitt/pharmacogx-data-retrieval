"""
PharmacoGx Data Prep Pipeline
V. Keith Hughitt
"""
import os

# output directory
output_dir = os.path.join(config['output_dir'], config['version'])
report_dir = os.path.join(config['report_dir'], config['version'])

# Generate lists of output filepath components
feature_datasets = []
pheno_datasets   = []

feature_types = []
phenotypes    = []

for dataset_id in config['datasets']:
    for feature_type in config['datasets'][dataset_id]['features']:
        feature_datasets.append(dataset_id)
        feature_types.append(feature_type)
    for phenotype in config['datasets'][dataset_id]['phenotypes']:
        pheno_datasets.append(dataset_id)
        phenotypes.append(phenotype)

rule all:
    input:
        expand(os.path.join(output_dir, '{dataset}/features/{feature_type}.feather'), 
               zip,
               dataset=feature_datasets, feature_type=feature_types),
        expand(os.path.join(output_dir, '{dataset}/features/{feature_type}_col_mdata.feather'), 
               zip,
               dataset=feature_datasets, feature_type=feature_types),
        expand(os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}.feather'), 
               zip,
               dataset=pheno_datasets, phenotype=phenotypes),
        expand(os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}_row_mdata.feather'), 
               zip,
               dataset=pheno_datasets, phenotype=phenotypes),
        expand(os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}_col_mdata.feather'), 
               zip,
               dataset=pheno_datasets, phenotype=phenotypes)

rule all_features:
    input:
        expand(os.path.join(output_dir, '{dataset}/features/{feature_type}.feather'), 
               zip,
               dataset=feature_datasets, feature_type=feature_types)

rule load_features:
    output:
        os.path.join(output_dir, '{dataset}/features/{feature_type}.feather'),
        os.path.join(output_dir, '{dataset}/features/{feature_type}_col_mdata.feather')
    script: 'src/load_features.R'

rule load_phenotypes:
    output:
        os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}.feather'),
        os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}_row_mdata.feather'),
        os.path.join(output_dir, '{dataset}/phenotypes/{phenotype}_col_mdata.feather')
    script: 'src/load_phenotypes.R'

rule download_metadata:
    output:
        os.path.join(config['raw_dir'], 'cell_annotation_table.csv'),
        os.path.join(config['raw_dir'], 'drug_annotation_table.csv')
    shell:
        """
        curl https://pharmacodb.pmgenomics.ca/download?cell_annotation=y -o {output[0]}
        curl https://pharmacodb.pmgenomics.ca/download?drug_annotation=y -o {output[1]}
        """
