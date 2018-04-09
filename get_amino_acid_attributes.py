import csv

def get_amino_acid_data():
    csv_file_name = "amino_acid_data.csv"

    Data = {}
    with open(csv_file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row["ID"]
            del row["ID"]
            Data[name] = row
    return Data
