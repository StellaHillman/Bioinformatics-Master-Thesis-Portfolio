#this script was used to number to convert the json output form the icn3d IgStrand numbering output to csv

import json
import csv
import os

def flatten_json(data, parent_key='', sep='_'):
    items = []
    for k, v in data.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_json(v, new_key, sep=sep).items())
        elif isinstance(v, list):
            for idx, sub_item in enumerate(v):
                items.extend(flatten_json(sub_item, f"{new_key}_{idx}", sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def convert_json_to_csv(json_filepath, csv_writer, protein_name):
    #read JSON file
    with open(json_filepath, 'r') as json_file:
        data = json.load(json_file)
    
    #flatten JSON structure
    flattened_data = flatten_json(data)
    
    #write rows to CSV with an additional column for the protein name
    for key, value in flattened_data.items():
        csv_writer.writerow([protein_name, key, value])

#input and output directories
input_directory = '/Users/stella/Desktop/Master_Projekt/IgStrand/FnIII-numbering/json_clean'
output_filepath = '/Users/stella/Desktop/Master_Projekt/IgStrand/FnIII-numbering/csv/merged-fnIII_output.csv'  # Change this path as needed

#open csv file to write all data into one file
with open(output_filepath, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    # Write header row with protein name, key, and value
    csv_writer.writerow(["Protein Name", "Key", "Value"])
    
    #iterate over each JSON file in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith('.json'):
            json_filepath = os.path.join(input_directory, filename)
            protein_name = os.path.splitext(filename)[0]  # Get protein name (filename without extension)
            print(f"Converting {json_filepath} to CSV with protein name {protein_name}")
            convert_json_to_csv(json_filepath, csv_writer, protein_name)

print(f"All files have been converted and saved in {output_filepath}")
