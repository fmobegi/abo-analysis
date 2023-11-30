#!/usr/bin/env python3
import os
import re
import sys
import xlsxwriter
import pandas as pd

class ABOReportParser:
    """
    A script to collate all ABO phenotype results from each sample into an 
    Excel wooksheet and a CSV for export to soft.
    """
    def __init__(self, input_dir):
        """
        Initialize the ABOReportParser.

        Args:
            input_dir (str): The input directory containing data files.
        """
        self.input_dir = input_dir
        self.results = []
        self.initialize_columns()

    def initialize_columns(self):
        """
        Define the column headers.
        """
        exon6 = ['Exon6_pos22'] * 10
        exon7_422 = ['Exon7_pos422'] * 10
        exon7_429 = ['Exon7_pos429'] * 10
        max_len = max(len(exon6), len(exon7_422), len(exon7_429))
        exon6 += [''] * (max_len - len(exon6))
        exon7_422 += [''] * (max_len - len(exon7_422))
        exon7_429 += [''] * (max_len - len(exon7_429))
        header_cols = ['', ''] + exon6 + exon7_422 + exon7_429 + ['', '', '']
        header_rows = ['Barcode', 'Sequencing_ID'] + ['#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type'] * 3 + ['Phenotype', 'Genotype', 'Reliability']
        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

    def parse_exon7(self, filename):
        """
        Open the file for reading
        """
        with open(filename, 'r', encoding="utf-8") as f:
            lines = f.readlines()
            pos1 = int(lines[2].split(':')[1].strip())
            # nuc1 = lines[3].split(':')[1].strip()
            # blood1 = lines[4].split(':')[1].strip()
            count1 = int(lines[6].split(':')[1].strip())
            mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]
            pos2 = int(lines[10].split(':')[1].strip())
            # nuc2 = lines[11].split(':')[1].strip()
            # blood2 = lines[12].split(':')[1].strip()
            count2 = int(lines[14].split(':')[1].strip())
            mat2, mis2, ins2, del2, a2, g2, c2, t2 = [int(x) for x in lines[16].split()]

        df = pd.DataFrame({
            'Exon': ['7', '7'],
            'Position': [pos1, pos2],
            '#Reads': [count1, count2],
            'Mat': [mat1, mat2],
            'Mis': [mis1, mis2],
            'Ins': [ins1, ins2],
            'Del': [del1, del2],
            'A': [a1, a2],
            'G': [g1, g2],
            'C': [c1, c2],
            'T': [t1, t2]
        })

        df['Position'] = df['Position'].apply(lambda x: x)
        
        def get_type(row):
            """Add the 'Type' column based on values in the [Ins Del A G C T] columns"""
            if row['Position'] == 422:
                if row['A'] >= 80:
                    return 'B'
                elif row['C'] >= 80:
                    return 'A or O'
                elif abs(row['A'] - row['C']) <= 20:
                    return '(A or O) and B'
                elif 20 < row['A'] < 80 and 20 < row['C'] < 80:
                    return '(A or O) and B'
                
            if row['Position'] == 429:
                if row['G'] >= 80:
                    return 'A or O'
                elif row['C'] >= 80:
                    return 'B'
                elif abs(row['G'] - row['C']) <= 20:
                    return '(A or O) and B'
                elif 20 < row['G'] < 80 and 20 < row['C'] < 80:
                    return '(A or O) and B'
            return ''

        df['Type'] = df.apply(get_type, axis=1)

        # Reorder the columns
        df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]

        return df

    def parse_exon6(self, filename):
        """Parse exon 6 and extract relevant data"""
        with open(filename, 'r', encoding="utf-8") as f:
            lines = f.readlines()

            # Extract the required lines from the first section
            pos1 = int(lines[2].split(':')[1].strip())
            count1 = int(lines[6].split(':')[1].strip())
            mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

        # Create a DataFrame from the extracted values
        df = pd.DataFrame({
            'Exon': ['6'],
            'Position': [pos1],
            '#Reads': [count1],
            'Mat': [mat1],
            'Mis': [mis1],
            'Ins': [ins1],
            'Del': [del1],
            'A': [a1],
            'G': [g1],
            'C': [c1],
            'T': [t1]
        })

        # Remove the string "Exon 6 pos" from the Position column
        df['Position'] = df['Position'].apply(lambda x: x)
        max_g = df['G'].max()
        max_del = df['Del'].max()
        if max_g >= 80 and max_g > max_del:
            df['Type'] = 'A or B'
        elif max_del >= 80 and max_del > max_g:
            df['Type'] = 'O'
        elif abs(max_g + max_del) >= 20:
            df['Type'] = 'O and (A or B)'
        elif 20 < df['Del'] < 80 and 20 < df['G']< 80:
            return 'O and (A or B)'
        else:
            df['Type'] = ''

        # Reorder the columns
        df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]

        return df

    def add_phenotype_genotype(self, df):
        """Find the column numbers with "Type" in the column name"""
        type_cols = [i for i, col in enumerate(df.columns) if 'Type' in col]

        # Assign the column numbers to Type1, Type2, and Type3
        Type1 = df.columns[type_cols[0]]
        Type2 = df.columns[type_cols[1]]
        Type3 = df.columns[type_cols[2]]

        # Initialize the Reliability column
        df['Reliability'] = ''

        for i in range(len(df)):
            if ((df.loc[i, Type1] == 'A or B').all() and (df.loc[i, Type2] == 'A or O').all() and (df.loc[i, Type3] == 'A or O')).any():
                df.at[i, 'Phenotype'] = 'A'
                df.at[i, 'Genotype'] = 'AA'
            else:
                df.at[i, 'Phenotype'] = ''
                df.at[i, 'Genotype'] = ''

        return df

    def assign_phenotype_genotype(self, df):
        """Assin the phenotype and genotype information """
        exon6 = df[('Exon6_pos22', 'Type')]
        exon7_422 = df[('Exon7_pos422', 'Type')]
        exon7_429 = df[('Exon7_pos429', 'Type')]

        nreads6 = df[('Exon6_pos22', '#Reads')]
        nreads_71 = df[('Exon7_pos422', '#Reads')]
        nreads_72 = df[('Exon7_pos429', '#Reads')]

        if (exon6 == 'A or B').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
            Phenotype = 'A'
            Genotype = 'AA'
            # Reliability = 'Enter-manually'

        elif (exon6 == 'A or B').all() and (exon7_422 == 'B').all() and (exon7_429 == 'B').all():
            Phenotype = 'B'
            Genotype = 'BB'
            # Reliability = 'Enter-manually'

        elif (exon6 == 'A or B').all() and (exon7_422 == '(A or O) and B').all() & (exon7_429 == '(A or O) and B').all():
            Phenotype = 'AB'
            Genotype = 'AB'
            # Reliability = 'Enter-manually'

        elif (exon6 == 'O').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
            Phenotype = 'O'
            Genotype = 'OO'
            # Reliability = 'Enter-manually'

        elif (exon6 == 'O and (A or B)').all() and (exon7_422 == 'A or O').all() and (exon7_429 == 'A or O').all():
            Phenotype = 'A'
            Genotype = 'AO'
            # Reliability = 'Enter-manually'

        elif (exon6 == 'O and (A or B)').all() and (exon7_422 == '(A or O) and B').all and (exon7_429 == '(A or O) and B').all():
            Phenotype = 'B'
            Genotype = 'BO'
            # Reliability = 'Enter-manually'
        else:
            Phenotype = 'Unknown'
            Genotype = 'Unknown'
            # Reliability = 'Enter-manually'

        if (nreads6 < 30).all() and (nreads_71 < 30).all() and (nreads_72 < 30).all():
            Reliability = 'Low #reads'
        else:
            Reliability = ''

        df[('', 'Phenotype')] = Phenotype
        df[('', 'Genotype')] = Genotype
        df[('', 'Reliability')] = Reliability

        return df

    def process_file(self, filename):
        sample_name, barcode = filename.split("_")
        exon6_dir = os.path.join(self.input_dir, filename, "exon6")
        exon7_dir = os.path.join(self.input_dir, filename, "exon7")

        # Check if both exon6 and exon7 directories exist
        if not (os.path.exists(exon6_dir) and os.path.exists(exon7_dir)):
            print(f"Skipping file {filename}. Missing exon6 or exon7 directory.")
            return

        exon6_file = os.path.join(exon6_dir, "ABOPhenotype.txt")
        exon7_file = os.path.join(exon7_dir, "ABOPhenotype.txt")

        # Check if both exon6 and exon7 ABOPhenotype.txt files exist
        if not (os.path.exists(exon6_file) and os.path.exists(exon7_file)):
            print(f"Skipping file {filename}. Missing ABOPhenotype.txt file in exon6 or exon7.")
            return

        # Define sample_df at the beginning of the method
        sample_df = pd.DataFrame({
            'Barcode': [barcode],
            'Sequencing_ID': [sample_name]
        })

        try:
            df_exon6 = self.parse_exon6(exon6_file)
        except Exception as e:
            print(f"Error processing exon6 for file {filename}: {str(e)}")
            df_exon6 = None

        try:
            df_exon7 = self.parse_exon7(exon7_file)
        except Exception as e:
            print(f"Error processing exon7 for file {filename}: {str(e)}")
            df_exon7 = None

        if df_exon6 is not None and df_exon7 is not None:
            df_exon7_pos422 = df_exon7.iloc[[0]].reset_index(drop=True)
            df_exon7_pos429 = df_exon7.iloc[[1]].reset_index(drop=True)
            merged_df = pd.concat([sample_df, df_exon6, df_exon7_pos422, df_exon7_pos429],
                                axis=1, join="inner").drop(['Exon', 'Position'], axis=1)
            merged_df['Barcode'] = merged_df['Barcode'].str.replace('barcode', '', case=False)
            merged_df['Barcode'] = pd.to_numeric(merged_df['Barcode'], errors='coerce')
            merged_df = self.add_phenotype_genotype(merged_df)
            merged_df.columns = self.columns
            merged_df = self.assign_phenotype_genotype(merged_df)
            self.results.append(merged_df)

    def process_files(self):
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                if match := re.match(
                    r"^(IMM|INGS|NGS)(-[0-9]+-[0-9]+)?_barcode\d+", filename
                ):
                    print("Processing file: " + filename)
                    # Extract barcode and sample_name from the filename
                    sample_name, barcode = filename.split("_")
                    self.process_file(filename)
                    print("Done adding Sample %s with barcode %s to merged data frame\n" % (sample_name, barcode))

    def merge_dataframes(self):
        final_df = pd.concat(self.results)
        final_df = final_df.sort_values(('', 'Barcode'))
        return final_df

    def save_results_to_file(self, final_df):
        # Write to text file
        final_df.to_csv('./ABO_result.txt', sep = '\t', index = False)
        writer = pd.ExcelWriter('./ABO_result.xlsx', engine = 'xlsxwriter')
        final_df.columns = final_df.columns.droplevel()
        final_df.to_excel(
            writer,
            sheet_name = 'ABO_Result',
            header = True,
            index = False,
            startrow = 1
        )
        
        workbook = writer.book
        worksheet = writer.sheets['ABO_Result']

        data_format = workbook.add_format({'bg_color': 'white', 'font_color': 'black', 'border': 1})
        header_format = workbook.add_format({'bold': True, 'fg_color': '#007399', 'border': 1, 'font_color': 'white'})
        red_bg_format = workbook.add_format({'bg_color': '#e2725b', 'font_color': 'black'})

        header_format.set_align('center')
        header_format.set_align('vcenter')

        num_rows, num_cols = final_df.shape

        try:
            worksheet.conditional_format(
                'A1:AI1000',
                {
                    'type': 'formula',
                    'criteria': '=$AI1="Low #reads"',
                    'format': red_bg_format,
                }
            )
        except Exception as e:
            print(f"An error occurred while applying conditional formatting: {str(e)}")

        for row in range(2, num_rows + 2):  # Add 2 to account for the header row
            for col in range(num_cols):
                cell_value = final_df.iat[row - 2, col]
                if not pd.isna(cell_value):
                    worksheet.write(row, col, cell_value, data_format)
        
        header_columns = ['Exon6_pos22', 'Exon7_pos422', 'Exon7_pos429']
        merge_ranges = [('C1:L1', header_columns[0]), ('M1:V1', header_columns[1]), ('W1:AF1', header_columns[2])]
        worksheet.merge_range('A1:B1', 'Sample', header_format)
        worksheet.merge_range('AG1:AI1', 'Result', header_format)


        
        for merge_range in merge_ranges:
            worksheet.merge_range(merge_range[0], merge_range[1], header_format)

        for col in range(num_cols):
            cell_value = final_df.columns[col]
            if not pd.isna(cell_value):
                worksheet.write(1, col, cell_value, header_format)

        writer.close()

        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""
        # self.df_for_lis_soft["ABO Geno Type1"] = final_df["Genotype"].str[0]
        # self.df_for_lis_soft["ABO Geno Type2"] = final_df["Genotype"].str[1]
        # Only process if not Unknown
        self.df_for_lis_soft["ABO Geno Type1"], self.df_for_lis_soft["ABO Geno Type2"] = (
            (final_df["Genotype"].str[0], final_df["Genotype"].str[1]) 
            if 'Genotype' in final_df.columns and (final_df['Genotype'] != 'Unknown').all()
            else ("", "")
        )

        self.df_for_lis_soft["ABO Pheno Type"] = final_df["Phenotype"]
        self.df_for_lis_soft["RH"] = ""
        self.df_for_lis_soft["Blood Type"] = final_df["Phenotype"]
        self.df_for_lis_soft["ABORH Comments"] = ""
        self.df_for_lis_soft.to_csv("./final_export.csv", index = False, encoding = "utf-8")

    def run(self):
        """
        Run the ABOReportParser.
        This method processes files, merges dataframes, and saves results.
        """
        self.process_files()
        final_df = self.merge_dataframes()
        print("Final Results:")
        print("-" * 216)
        print(final_df.to_string(index=False))
        print("-" * 216)
        self.save_results_to_file(final_df)

if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) != 2:
        print("\nUsage: python ABOReportParser.py <input_directory>\n")
        sys.exit(1)
    input_directory = sys.argv[1]
    parser = ABOReportParser(input_directory)
    parser.run()
    print("""All done!\n""")

sys.exit(0)
