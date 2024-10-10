import os
import re

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

"""
Dear dev, 

here you'll find helper functions to load the spectra data and 
the associated metadata from circular dichroism data files

For each data loading function, the following values should be returned:

    wavelength (vector), spectra (matrix), spectra names (vector), 
    high tension voltage (matrix) 
 
"""

def file_ends_with_pattern(file_path):
    pattern = r"\.d\d+$"  # Regular expression pattern to match '.d' followed by one or more digits

    # Use re.search() to check if the pattern is found at the end of the file path
    if re.search(pattern, file_path):
        return True
    else:
        return False


def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False


def is_ht_curve(string):
    string = string.lower()

    c1 = "ht" in string
    c2 = "vt" in string
    c3 = "dynode" in string
    c4 = "volt" in string
    c5 = "dynode" in string

    return any([c1, c2, c3, c4, c5])


def column_name_is_important(string):
    string = string.lower()

    c1 = "lambda" in string
    c2 = "wavelength" in string
    c3 = "cd" in string
    c4 = "volts" in string
    c5 = "temperature" in string

    return any([c1, c2, c3, c4, c5])


def are_all_strings_numeric(lst):
    for item in lst:
        if not all(char.isdigit() or char in [".", "-"] for char in item):
            return False
    return True


# Guess if a string contains ',' separated values or ';' separated values
# Return either ',' or ';'
def find_delimiter_character(string):
    comma_count = string.count(",")
    semicolon_count = string.count(";")

    if semicolon_count > comma_count:

        return ";"

    else:

        return ","


def file_is_chirakit_txt_with_header(file):
    with open(file, encoding="latin-1") as f:

        ls = f.read().splitlines()

        # Iterate over all the lines except the last 10. We assume that we need at least 9 data points...
        for i, l in enumerate(ls[:-10]):

            if ("Wavelength_(nm)" in l or "Wavelength (nm)" in l) and any(
                ["CD_/_" in x for x in l.split()]
            ):

                next_line_elements = ls[i + 1].split()
                next_line_elements = [item for item in next_line_elements if item != ""]

                if (
                    are_all_strings_numeric(next_line_elements)
                    and len(next_line_elements) > 1
                ):
                    return True

    return False


def file_is_jasco_single_sample_csv(file):
    with open(file, encoding="latin-1") as f:

        ls = f.read().splitlines()

        for l in ls[:20]:

            l_semicolon_split = l.split(";")
            l_comma_split = l.split(",")

            condition1 = (
                l_semicolon_split[0].lower() == "origin"
                and l_semicolon_split[1].lower() == "jasco"
            )
            condition2 = (
                l_comma_split[0].lower() == "origin"
                and l_comma_split[1].lower() == "jasco"
            )

            if condition1 or condition2:
                return True

    return False

def read_jasco_single_sample_csv(file):
    with open(file, encoding="latin-1") as f:

        ls = f.read().splitlines()

        split_str = find_delimiter_character("".join(ls[:20]))

        # Find ydata description

        y_data_init_idx = [
            i for i, l in enumerate(ls[:20]) if l.split(split_str)[0] == "XUNITS"
        ][0] + 1
        y_data_end_idx = [
            i for i, l in enumerate(ls[:20]) if l.split(split_str)[0] == "FIRSTX"
        ][0]

        y_data_names = [l.split(split_str)[1] for l in ls[y_data_init_idx:y_data_end_idx]]

        n_cols = len(y_data_names) + 1

        # Initialize a list to store CD data for each column
        cd_data = [[] for _ in range(n_cols - 1)]
        # Initialize an empty list to store wavelengths
        wavelength = []

        # Iterate through the lines containing CD data (starting from row i+2)
        for i2, l in enumerate(ls[21:]):
            # Split the line using comma as the delimiter and remove any empty elements
            row = l.split(split_str)
            row = [x for x in row if x]

            # replace with dots, if needed
            if split_str == ";":
                row = [x.replace(",", ".") for x in row]

            # Check if all elements in the row are float values
            if not all([is_float(x) for x in row]):
                break

            # Check if the row has the correct number of columns (nCols)
            if len(row) != n_cols:
                break

            # Append the wavelength value to the list
            wavelength.append(float(row[0]))
            # Append the CD data for each column to the corresponding list in cdData
            for i3, r in enumerate(row[1:]):
                cd_data[i3].append(float(r))

        # Convert the lists in cdData to numpy arrays
        cd_data = [np.array(x) for x in cd_data]

        # Create a DataFrame from the CD data arrays
        df = pd.DataFrame.from_records(cd_data)

        # Transpose the DataFrame to have spectra as rows
        signal = df.transpose()

        # Assign appropriate column names to the signal DataFrame
        signal.columns = y_data_names

        # Add a new 'wavelength' column to the signal DataFrame
        signal["wavelength"] = wavelength

        cd_signal_name = [x for x in y_data_names if "cd" in x.lower()][0]

        wavelength = np.array(wavelength)

        spectra = np.array(signal[[cd_signal_name]]).astype("float")

        # Create an array with column names for spectra
        spectra_names = np.array(os.path.basename(file))

        try:

            ht_signal_name = [x for x in y_data_names if is_ht_curve(x)][0]
            signal_ht = np.array(signal[[ht_signal_name]]).astype("float")

        except:

            signal_ht = np.empty_like(spectra)
            signal_ht[:] = np.nan  # Fill the array with NaN values

        return wavelength, spectra, spectra_names, signal_ht

def read_jasco_single_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="latin-1") as f:

        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find delimiter
        split_str = find_delimiter_character("".join(ls[:10]))

        # Avoid consecutive ';'' or ','
        ls = [re.sub(";+", ";", l) for l in ls]
        ls = [re.sub(",+", ",", l) for l in ls]

        # Find the row index where the numeric data starts (CD spectra)
        for l in ls:

            l2 = l.split(split_str)
            l2 = [x for x in l2 if x != ""]

            # Check if the element has both key and value (length should be 2)
            if len(l2) == 2:
                # Add the key-value pair to the metadata dictionary
                metadata[l2[0]] = l2[1]

    # Return the extracted metadata dictionary
    return metadata

def detect_file_type(file):
    # Define a dictionary to map file extensions to their corresponding types
    file_extensions = {
        ".gen": "GenCDfile",  # If the file has the extension '.gen', it is of type 'GenCDfile'
        ".pcd": "PCCDBFile",  # If the file has the extension '.pcd', it is of type 'PCCDBFile'
        ".dat": "DatFile",  # If the file has the extension '.dat', it is of type 'DatFile'
    }

    # Check if the file matches any of the specified file extensions
    for ext, value in file_extensions.items():
        if file.endswith(ext):
            return value

    # Check if the file ends with a pattern '.d' followed by one or more digits (using the file_ends_with_pattern function)
    if file_ends_with_pattern(file):
        return "d0xFile"  # If the pattern is found, the file is of type 'd0x'

    if file_is_jasco_single_sample_csv(file):
        return "jasco_simple"

    # If the file doesn't match any of the previous cases, open the file and read its lines
    with open(file, encoding="latin-1") as f:
        ls = f.read().splitlines()
        # Iterate through each line and check if it contains both "Wavelength" and "Temperature"
        for l in ls:
            if "Wavelength" in l and "Temperature" in l:
                return "ChirascanFileTemperatureRamp"  # If both keywords are found, the file is of type 'ChirascanFileTemperatureRamp'

    if file_is_chirakit_txt_with_header(file):
        return "chirakit_txt_with_header"

    # If the file doesn't match any of the previous cases, open the file and read its lines
    with open(file, encoding="latin-1") as f:
        ls = f.read().splitlines()
        # Iterate through each line and check if it contains both "Wavelength" and "Temperature"
        for i, l in enumerate(ls):

            if "Wavelength" in l and "CircularDichroism" in ls[i - 1]:
                return "ChirascanFile"  # If both keywords are found, the file is of type 'ChirascanFileTemperatureRamp'

    # If none of the previous conditions are met, the file type is 'ChirascanFile'
    return "plain_csv"


def read_custom_csv(file):
    with open(file, encoding="latin-1") as f:

        ls = f.read().splitlines()[:80]

        # find if the delimitr could be a comma or a semicolon
        split_str = find_delimiter_character("".join(ls))

    # Try comma or semicolon delimiters.
    df = pd.read_csv(file, comment="#", encoding="latin1", delimiter=split_str)

    # Check if all columns are numeric
    all_numeric = all(is_numeric_dtype(df[col]) for col in df.columns)

    # try changing the decimal to ','
    if not all_numeric:
        df = pd.read_csv(
            file, comment="#", encoding="latin1", delimiter=split_str, decimal=","
        )

    # If we have only one column, change separator delimiter to spaces
    if df.shape[1] == 1:

        df = pd.read_csv(file, comment="#", delimiter=r"\s+", encoding="latin1")

        all_numeric = all(df[col].dtype == np.number for col in df.columns)

        # try changing the decimal to ','
        if not all_numeric:
            df = pd.read_csv(
                file, comment="#", encoding="latin1", delimiter=r"\s+", decimal=","
            )

    if df.shape[1] == 3:

        df_wide = df.pivot_table(
            index=df.columns[0], columns=df.columns[1], values=df.columns[2]
        ).reset_index()

        wavelength = np.array(df_wide.iloc[:, 0])
        spectra = np.array(df_wide.iloc[:, 1:], dtype="float")
        spectra_names = np.array(df_wide.columns[1:], dtype="str")

    elif df.shape[1] == 2:

        file_name = os.path.basename(file)

        wavelength = np.array(df.iloc[:, 0])

        # We have one CD spectrum
        if len(np.unique(wavelength)) == len(wavelength):

            spectra = np.array(
                df.iloc[:, 1:], dtype="float"
            )  # don't use np.array(df.iloc[:,1]) because we need to keep the same dimensions ...
            spectra_names = np.array([file_name], dtype="str")

        # We have the CD signal measured at one (or more) wavelength(s) as a
        # function of an unkown parameter
        else:

            spectra = []
            unique_wls = np.unique(wavelength)

            for wl in unique_wls:
                
                temp_df = df.iloc[np.where(wavelength == wl)]
                spectra.append(temp_df.values[:, 1])

            spectra = np.array(spectra, dtype="float")
            spectra_names = np.array(
                [str(i) + " " + file_name for i in range(spectra.shape[1])], dtype="str"
            )
            wavelength = unique_wls
    else:

        wavelength = np.array(df.iloc[:, 0])
        spectra = np.array(df.iloc[:, 1:], dtype="float")
        spectra_names = np.array(df.columns[1:], dtype="str")

    signal_ht = np.empty_like(spectra)
    signal_ht[:] = np.nan  # Fill the array with NaN values

    return wavelength, spectra, spectra_names, signal_ht

def read_chirakit_txt_data(file):
    df = pd.read_csv(file, comment="#", delimiter=r"\s+", encoding="latin1")

    # Find if there is voltage signal
    contains_voltage = any("HT_/_" in s for s in df.columns)

    wavelength = np.array(df.iloc[:, 0])

    last_cd_signal_column_id = (
        int((len(df.columns) - 1) - (contains_voltage * (len(df.columns) - 1) / 2)) + 1
    )

    spectra = np.array(df.iloc[:, 1:last_cd_signal_column_id])

    if contains_voltage:

        signal_ht = np.array(df.iloc[:, last_cd_signal_column_id:])

    else:

        signal_ht = np.empty_like(spectra)
        signal_ht[:] = np.nan  # Fill the array with NaN values

    spectra_names = df.columns[1:last_cd_signal_column_id]
    spectra_names = [s.replace("CD_/_", "") for s in spectra_names]
    spectra_names = np.array(spectra_names)

    return wavelength, spectra, spectra_names, signal_ht

def read_chirakit_txt_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="utf-8") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Extract lines that start with '#' and remove the leading ';'
        info = [l[1:] for l in ls if l.startswith("#")]

        # Extract metadata before the "Comment" section
        for l in info:
            # Split each line into key and value using fixed column positions
            key = " ".join(l[:47].split())
            value = " ".join(l[47:].split())

            # Remove the last two characters which should be ' :'
            key = key[:-2]

            # Add the key-value pair to the metadata dictionary
            metadata[key] = value

    # Return the extracted metadata dictionary
    return metadata

def extract_d0x_data(file):
    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find where the numeric data starts
        for i, l in enumerate(ls):

            l2 = l.lower()
            c1 = "lambda" in l2
            c2 = "wavelength" in l2
            c3 = "cd" in l2

            if (c1 or c2) and c3:
                break

        # Assign column names
        col_names = ls[i][1:].split()

        selected_colnames = [x for x in col_names if column_name_is_important(x)]

        dataset_lines = ls[(i + 1) :]
        df = pd.DataFrame([line.split() for line in dataset_lines], columns=col_names)

        df = df[selected_colnames]

        df[selected_colnames] = df[selected_colnames].astype(float)

        return df

def read_d0x_file_data(file):
    df = extract_d0x_data(file)

    # select the wavelength (1st column) and spectral (2nd and 3rd columns) data
    wavelength = np.array(df.iloc[:, 0])
    spectra = np.array(df.iloc[:, 1])
    signal_ht = np.array(df.iloc[:, 2])
    spectra_names = np.array(["CD mdeg"])

    # Convert 1D array to column vector
    signal_ht = signal_ht[:, np.newaxis]
    spectra = spectra[:, np.newaxis]

    # Return the wavelength array, spectra data, and spectrum names
    return wavelength, spectra, spectra_names, signal_ht


def read_temperature_d0x(file):
    df = extract_d0x_data(file)
    temperature = np.array(df.iloc[1, 3])

    return temperature


def read_d0x_file_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Extract lines that start with ';' and remove the leading ';'
        info = [l[1:] for l in ls if l.startswith(";")]

        # Find the index where the "Comment" section starts
        for i, l in enumerate(info):
            if "Comment" in l:
                break

        # Extract metadata before the "Comment" section
        for l in info[:i]:
            # Split each line into key and value using fixed column positions
            key = " ".join(l[:30].split())
            value = " ".join(l[30:].split())

            # Add the key-value pair to the metadata dictionary
            metadata[key] = value

        # Extract metadata after the "Comment" section
        for l in info[i:]:
            # Split each line by ':' to get key and value
            splitted = l.split(":")
            splitted = [item for item in splitted if item != ""]

            # If the line has exactly two parts, consider it as a valid key-value pair
            if len(splitted) == 2:
                key = " ".join(splitted[0].split())
                value = " ".join(splitted[1].split())

                # Add the key-value pair to the metadata dictionary
                metadata[key] = value

    # Return the extracted metadata dictionary
    return metadata


def read_dat_file_data(file):
    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the indices where the data starts and ends
        start_indices = [
            index + 2 for index, line in enumerate(ls) if line.startswith("$MDCDATA")
        ]
        end_index = [
            index for index, line in enumerate(ls) if line.startswith("$ENDDATA")
        ]
        end_indices = [s - 3 for s in start_indices[1:]] + end_index

        # Create an empty list to store DataFrames
        dfs = []
        dfs_ht = []
        # Loop through the identified start and end indices to extract and process individual spectra
        for s, e in zip(start_indices, end_indices):

            # Get the scan name and column names for the DataFrame
            scan_name = ls[s - 3].split(":")[1]
            col_names = ["wavelength"] + [
                scan_name + " " + x for x in ls[s - 1].split()[1:]
            ]

            # Extract the lines corresponding to the current spectrum
            dataset_lines = ls[s:e]

            # Create a DataFrame from the extracted data and convert to numeric
            df = pd.DataFrame(
                [line.split() for line in dataset_lines], columns=col_names
            )
            df[col_names] = df[col_names].astype(float)

            # Remove unwanted columns containing 'Jacket_Temp' or 'Error'
            idx_to_remove = [x for x in col_names if "Jacket_Temp" in x or "Error" in x]
            df = df.drop(idx_to_remove, axis=1)

            # Sort the DataFrame based on the 'wavelength' column
            df.sort_values("wavelength", inplace=True)

            # Create the HT voltage dataframe
            # Extract the voltage curve
            ht_curve = [s for s in df.columns if is_ht_curve(s)]
            if len(ht_curve) > 0:

                column_to_repeat = ht_curve[0]

                repeated_column = [
                    df[column_to_repeat].values for _ in range(len(df.columns) - 2)
                ]
                repeated_column = np.array(repeated_column).flatten()
                repeated_column_matrix = repeated_column.reshape(
                    -1, len(df.columns) - 2
                )

                df = df.drop(column_to_repeat, axis=1)
                df_signal_ht = pd.DataFrame(
                    repeated_column_matrix, columns=df.columns[1:]
                )
                df_signal_ht.insert(0, "wavelength", df["wavelength"].values)

            else:

                df_signal_ht = pd.DataFrame(
                    np.nan, index=range(df.shape[0]), columns=range(df.shape[1])
                )

            # Append the processed DataFrame to the list of DataFrames
            dfs.append(df)
            dfs_ht.append(df_signal_ht)

        # Case 1 - only one scan
        if len(dfs) == 1:
            merged = dfs[0]
            merged_ht = dfs_ht[0]
        # Case 2 - two or more spectra
        else:
            # Merge the 1st and second spectra data
            if len(dfs) > 1:
                merged = pd.merge_asof(
                    dfs[0].dropna(),
                    dfs[1].dropna(),
                    on="wavelength",
                    direction="nearest",
                    allow_exact_matches=True,
                )
                merged_ht = pd.merge_asof(
                    dfs_ht[0].dropna(),
                    dfs_ht[1].dropna(),
                    on="wavelength",
                    direction="nearest",
                    allow_exact_matches=True,
                )

            # Merge the remaining spectral data
            if len(dfs) > 2:
                for df, df_ht in zip(dfs[2:], dfs_ht[2:]):
                    merged = pd.merge_asof(
                        merged,
                        df.dropna(),
                        on="wavelength",
                        direction="nearest",
                        allow_exact_matches=True,
                    )
                    merged_ht = pd.merge_asof(
                        merged_ht,
                        df_ht.dropna(),
                        on="wavelength",
                        direction="nearest",
                        allow_exact_matches=True,
                    )

        # Get the list of spectrum names and convert the merged DataFrame to a numpy array
        spectra_names = merged.columns.tolist()
        merged = np.array(merged, dtype="float")
        merged_ht = np.array(merged_ht, dtype="float")

        spectra = merged[:, 1:]
        signal_ht = merged_ht[:, 1:]
        wavelength = merged[:, 0]
        spectra_names = spectra_names[1:]

        # Return the wavelength array, spectra data, and spectra names
        return wavelength, spectra, spectra_names, signal_ht


def read_dat_file_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the index where metadata ends (before the '$MDCDATA' line)
        end_index = [
            index - 3 for index, line in enumerate(ls) if line.startswith("$MDCDATA")
        ][0]

        # Iterate through the lines containing metadata (before the '$MDCDATA' line)
        for l in ls[:end_index]:
            # Split the line using comma as the delimiter and remove any empty elements
            l2 = l.split(",")
            l2 = [x for x in l2 if x != ""]

            # Iterate through the elements in the split line
            for element in l2:
                # Split each element using colon as the delimiter and remove any leading/trailing whitespaces and empty elements
                splitted = element.split(":", 1)
                splitted = [s.strip() for s in splitted]
                splitted = [x for x in splitted if x != ""]

                # Check if the element has both key and value (len should be 2)
                if len(splitted) == 2:
                    # Add the key-value pair to the metadata dictionary
                    metadata[splitted[0]] = splitted[1]

    # Return the extracted metadata dictionary
    return metadata


def read_pccdb_file_data(file):
    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the row index where the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            if "DATA" in l and "Wavelength" in l:
                data_begin = i
                break

        # Extract the spectrum names from the 'DATA' row
        spectra_names = re.split(r"\d+\.", ls[data_begin])[1:]
        spectra_names = [s.strip() for s in spectra_names]
        spectra_names = [s[:-1] if s.endswith(")") else s for s in spectra_names]
        spectra_names = [s[:-1] if s.endswith(".") else s for s in spectra_names]

        # Extract the column values for each spectrum and exclude rows with 'CALIBRATION'
        column_values = [
            l.split()
            for l in ls[data_begin:]
            if len(l.split()) == len(spectra_names) and "CALIBRATION" not in l
        ]

        # Convert the list of column values to a numpy array (spectrum data)
        spectra = np.array(column_values, dtype="float")

        # Extract the wavelength array from the first column of the spectrum data
        wavelength = spectra[:, 0].flatten()

        spectra_names = spectra_names[1:]  # (excluding the first element)
        spectra = spectra[:, 1:]  # (excluding the first column)

        # Extract the voltage curve
        ht_index = [i for i, s in enumerate(spectra_names) if is_ht_curve(s)]
        if len(ht_index) > 0:

            signal_ht = spectra[:, ht_index[0]]

            spectra = np.delete(spectra, ht_index[0], axis=1)
            del spectra_names[ht_index[0]]

            # Repeat the 1D vector along columns to create a 2D array
            signal_ht = np.repeat(signal_ht, spectra.shape[1]).reshape(
                -1, spectra.shape[1]
            )

        else:
            signal_ht = np.empty_like(spectra)
            signal_ht[:] = np.nan  # Fill the array with NaN values

    # Return the wavelength vector, spectra array , and spectrum names vector
    return wavelength, spectra, spectra_names, signal_ht


def read_pccdb_file_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Iterate until  the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            if "DATA" in l and "Wavelength" in l:
                break

            # Extract the key and value from the line
            key = l[:60].strip()
            value = l[60:].strip()

            # Check if the value is not empty
            if value:
                # Add the key-value pair to the metadata dictionary
                metadata[key] = value

    # Return the extracted metadata dictionary
    return metadata

def read_chirascan_file_meta_data(cd_file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(cd_file) as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the row index where the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            # Check if the line contains "Wavelength"
            if "Wavelength" in l:
                # Assuming metadata starts two lines below the line containing "Wavelength"
                # Get the third row (i+2) and split it by comma to count the number of columns
                row2 = ls[i + 2].split(",")
                # Check if the first element of row2 is numeric, which indicates the start of data
                if row2[0].isnumeric():
                    break

        # Iterate through lines before the start of data (index i+2)
        for l in ls[: i + 2]:
            # Split the line using comma as the delimiter and remove any empty elements
            l2 = l.split(",")
            l2 = [x for x in l2 if x != ""]

            # Iterate through the elements in the split line
            for element in l2:
                # Split each element using colon as the delimiter and remove any leading/trailing whitespaces and empty elements
                splitted = element.split(":", 1)
                splitted = [s.strip() for s in splitted]
                splitted = [x for x in splitted if x != ""]

                # Check if the element has both key and value (length should be 2)
                if len(splitted) == 2:
                    # Add the key-value pair to the metadata dictionary
                    metadata[splitted[0]] = splitted[1]

    # Return the extracted metadata dictionary
    return metadata


def read_chirascan_file_data(cd_file):
    # Open the file for reading
    with open(cd_file) as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the row index where the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            # Check if the line contains "Wavelength"
            if "wavelength" in l.lower():
                # Assuming metadata starts two lines below the line containing "Wavelength"
                # Get the third row (i+2) and split it by comma to count the number of columns
                row2 = ls[i + 2].split(",")
                # Check if the first element of row2 is numeric, which indicates the start of data
                if row2[0].isnumeric():
                    # n_cols indicates the number of columns in the data
                    n_cols = len([x for x in row2 if x])
                    break

        # Initialize a list to store CD data for each column
        cd_data = [[] for _ in range(n_cols - 1)]
        # Initialize an empty list to store wavelengths
        wavelength = []

        # Iterate through the lines containing CD data (starting from row i+2)
        for i2, l in enumerate(ls[i + 2 :]):
            # Split the line using comma as the delimiter and remove any empty elements
            row = l.split(",")
            row = [x for x in row if x]

            # Check if all elements in the row are float values
            if not all([is_float(x) for x in row]):
                break

            # Check if the row has the correct number of columns (n_cols)
            if len(row) != n_cols:
                break

            # Append the wavelength value to the list
            wavelength.append(float(row[0]))
            # Append the CD data for each column to the corresponding list in cd_data
            for i3, r in enumerate(row[1:]):
                cd_data[i3].append(float(r))

        # Convert the lists in cd_data to numpy arrays
        cd_data = [np.array(x) for x in cd_data]

        # Create a DataFrame from the CD data arrays
        df = pd.DataFrame.from_records(cd_data)

        # Transpose the DataFrame to have spectra as rows
        signal = df.transpose()

        # Convert the list of wavelengths to a numpy array
        wavelength = np.array(wavelength).astype(float)

        # Assign appropriate column names to the signal DataFrame
        signal.columns = ["C" + str(x) for x in range(signal.shape[1])]

        # Add a new 'wavelength' column to the signal DataFrame
        signal["wavelength"] = wavelength

        # Group the signal DataFrame by 'wavelength' and calculate the mean for each group
        aggs = signal.groupby("wavelength")[signal.columns.values].mean()

        # Extract the mean wavelengths and spectra as numpy arrays
        wavelength = np.array(aggs.iloc[:, -1]).astype("float").flatten()
        spectra = np.array(aggs.iloc[:, :-1]).astype("float")

        # Create an array with column names for spectra
        spectra_names = np.array(["Col " + str(x + 1) for x in range(spectra.shape[1])])

        # Create an empty array for the HT values
        signal_ht = np.empty_like(spectra)
        signal_ht[:] = np.nan  # Fill the array with NaN values

    # Return the wavelength array, spectra data, and spectrum names
    return wavelength, spectra, spectra_names, signal_ht


def read_chirascan_file_data_thermal_ramp(cd_file):
    # Open the file for reading
    with open(cd_file) as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the row index where the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            # Check if the line contains "Wavelength"
            if "Wavelength" in l:
                # Assuming metadata starts two lines below the line containing "Wavelength"
                # Get the third row (i+2) and split it by comma to count the number of columns
                row2 = ls[i + 2].split(",")
                # Check if the first element of row2 is numeric, which indicates the start of data
                if row2[0].isnumeric():
                    # n_cols indicates the number of columns in the data
                    n_cols = len([x for x in row2 if x])
                    break

        # Initialize a list to store CD data for each column
        cd_data = [[] for _ in range(n_cols - 1)]
        # Initialize an empty list to store wavelengths
        wavelength = []

        # Extract temperature values from the line after the "Wavelength" line
        temperature = ls[i + 1].split(",")[1:]
        temperature = [t for t in temperature if t]
        # Convert temperature values to a numpy array
        temperature = np.array(temperature).astype(float)

        # Iterate through the lines containing CD data (starting from row i+2)
        for i2, l in enumerate(ls[i + 2 :]):
            # Split the line using comma as the delimiter and remove any empty elements
            row = l.split(",")
            row = [x for x in row if x]

            # Check if all elements in the row are float values
            if not all([is_float(x) for x in row]):
                break

            # Check if the row has the correct number of columns (n_cols)
            if len(row) != n_cols:
                break

            # Append the wavelength value to the list
            wavelength.append(float(row[0]))
            # Append the CD data for each column to the corresponding list in cd_data
            for i3, r in enumerate(row[1:]):
                cd_data[i3].append(float(r))

        # Convert the lists in cd_data to numpy arrays
        cd_data = [np.array(x) for x in cd_data]

        # Create a DataFrame from the CD data arrays
        df = pd.DataFrame.from_records(cd_data)
        # Transpose the DataFrame to have spectra as rows
        signal = df.transpose()

        # Convert the list of wavelengths to a numpy array
        wavelength = np.array(wavelength).astype(float)

        # Assign appropriate column names to the signal DataFrame
        signal.columns = ["C" + str(x) for x in range(signal.shape[1])]

        # Add a new 'wavelength' column to the signal DataFrame
        signal["wavelength"] = wavelength

        # Group the signal DataFrame by 'wavelength' and calculate the mean for each group
        aggs = signal.groupby("wavelength")[signal.columns.values].mean()

        # Extract the mean wavelengths and spectra as numpy arrays
        wavelength = np.array(aggs.iloc[:, -1]).astype("float").flatten()
        spectra = np.array(aggs.iloc[:, :-1]).astype("float")

        # Create an empty array for the HT values
        signal_ht = np.empty_like(spectra)
        signal_ht[:] = np.nan  # Fill the array with NaN values

    # Return the wavelength array, spectra data, and temperature values
    return wavelength, spectra, temperature, signal_ht


def read_gen_file_meta_data(file):
    # Initialize an empty dictionary to store metadata
    metadata = {}

    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Iterate through the lines in the file
        for l in ls:
            # Split the line using tab as the delimiter
            l2 = l.split("\t")

            # If the line has more than 3 elements, exit the loop
            # (assuming the metadata section ends at this point)
            if len(l2) > 3:
                break

            # Remove leading/trailing whitespaces from each element in the split line
            l3 = [s.strip() for s in l2]

            # Create a key-value pair for metadata (assuming the first element is the key, and the rest are the value)
            metadata[l3[0]] = " ".join(l3[1:])

    # Return the extracted metadata dictionary
    return metadata


def read_gen_file_data(file):
    # Open the file for reading
    with open(file, encoding="latin-1") as f:
        # Read all lines and split them into a list of lines
        ls = f.read().splitlines()

        # Find the row index where the numeric data starts (CD spectra)
        for i, l in enumerate(ls):
            # Split the line using tab as the delimiter
            l2 = l.split("\t")

            # If the line has more than 3 elements, exit the loop
            # (assuming the numeric data section starts at this point)
            if len(l2) > 3:
                break

        # Create a numpy array of column values for the numeric data (skipping the metadata section)
        column_values = np.array([x.split() for x in ls[i:]], dtype="float")
        # Define spectra names based on the expected structure
        spectra_names = np.array(
            [
                "Final Processed Spectrum",
                "HT values (from raw spectrum 1)",
                "Final Processed Spectrum",
                "Blank",
                "Blank",
                "Blank",
            ]
        )

        # Extract the wavelength column from the column values
        wavelength = column_values[:, 0].flatten()
        # Extract the spectra data (excluding the first column, which is wavelength)
        spectra = column_values[:, 1:]

        total_points = spectra.shape[0]
        idx_to_keep = []

        # Filter duplicated spectra without data and spectra with too many zeros or NaNs
        for i in range(spectra.shape[1]):
            duplicated = False
            count_nas = np.sum(np.isnan(spectra[:, i]))
            count_zeros = np.sum(np.where(spectra[:, i] == 0))

            # Avoid the iteration for the last column
            if i != spectra.shape[1]:
                # Check for duplicates among the remaining columns
                for ii in range(i + 1, spectra.shape[1]):
                    if np.array_equal(spectra[:, i], spectra[:, ii]):
                        duplicated = True
                        break

            # If the column is not duplicated and has an acceptable amount of data, keep it
            if (count_zeros + count_nas < (total_points / 4)) and not duplicated:
                idx_to_keep.append(i)

    spectra = spectra[:, idx_to_keep]
    spectra_names = spectra_names[idx_to_keep]

    # Extract the voltage curve
    ht_index = [i for i, s in enumerate(spectra_names) if is_ht_curve(s)]
    if len(ht_index) > 0:

        signal_ht = spectra[:, ht_index[0]]

        spectra = np.delete(spectra, ht_index[0], axis=1)
        spectra_names = np.delete(spectra_names, ht_index[0])

        # Repeat the 1D vector along columns to create a 2D array
        signal_ht = np.repeat(signal_ht, spectra.shape[1]).reshape(-1, spectra.shape[1])

    else:

        signal_ht = np.empty_like(spectra)
        signal_ht[:] = np.nan  # Fill the array with NaN values

    # Return the wavelength array, filtered spectra data, and the corresponding spectra names
    return wavelength, spectra, spectra_names, signal_ht

def read_unfolding_file_temperature(file):

    # Open the file for reading
    with open(file, encoding="utf-8") as f:
        # Read the first n lines and split them into a list of lines
        ls = f.read().splitlines()[:30]

        # Extract lines that start with '#' and remove the leading ';'
        info = [l[1:] for l in ls if l.startswith("#")]

        # Extract metadata before the "Comment" section
        for l in info:
            # Split each line into key and value using fixed column positions

            if (':' in l):

                key   = l.split(':')[0]
                value = l.split(':')[1]

                if 'temperature' in key.lower():

                    return float(re.sub(r'[^0-9]', '', value))

    return 0

def are_headers_numeric(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path, nrows=0,comment='#')  # Read only the headers

    # Extract the headers
    headers = df.columns[1:]

    # Check if each header can be interpreted as numeric
    for header in headers:
        try:
            float(header)
        except ValueError:
            return False
    return True

def read_unfolding_data_monomer(file_path):
    """
    Load a CSV file with unfolding data.

    The columns are ordered as follows:
    - The first column contains the temperature / chemical agent concentration.
    - Subsequent columns contain the signal of each curve.

    Args:
    - file_path (str): Path to the CSV file.

    Returns:
    - pd.DataFrame: A pandas DataFrame containing the temperature and signal data.
    """

    df = pd.read_csv(file_path, comment="#")

    signals_all       = []
    mment_factor_all  = []

    wavelength_useful = np.array(df.columns[1:])

    df.sort_values(by=df.columns[0], inplace=True)

    for i in range(1, len(df.columns)):

        df_temp = np.array(df.iloc[:, [0, i]])
        df_temp = df_temp[~np.isnan(df_temp).any(axis=1)]

        mment_factor_all.append(df_temp[:, 0])
        signals_all.append(df_temp[:, 1])

    mment_factor_all = np.concatenate(mment_factor_all)
    #signals          = np.array(signals_all).reshape(len(df.columns), -1)

    signal_matrix = np.full((len(wavelength_useful), len(mment_factor_all)), np.nan)

    n = 0

    for row in range(len(wavelength_useful)):

        m = n + len(signals_all[row])

        signal_matrix[row,n:m] = np.array(signals_all[row])

        n += m

    return mment_factor_all, signal_matrix, wavelength_useful

def read_unfolding_data_oligomer(file_path):
    """
    Load a CSV file with unfolding data.

    The columns are ordered as follows:
    - The first column contains the temperature / chemical agent concentration.
    - Subsequent columns contain the signal of each curve.
    - The header for each curve (except the temperature column) corresponds to the n-mer concentration in micromolar units.

    Args:
    - file_path (str): Path to the CSV file.

    Returns:
    - pd.DataFrame: A pandas DataFrame containing the temperature and signal data.
    """

    df = pd.read_csv(file_path, comment="#")

    concentrations = [re.sub(r"[^\d.]", "", col) for col in df.columns[1:]]

    mment_factor_all = []
    prot_conc_all    = []
    signals_all      = []

    wavelength_useful = np.array([""])

    df.sort_values(by=df.columns[0], inplace=True)

    for i in range(1, len(df.columns)):

        df_temp = np.array(df.iloc[:, [0, i]])
        df_temp = df_temp[~np.isnan(df_temp).any(axis=1)]

        mment_factor_all.append(df_temp[:, 0])
        signals_all.append(df_temp[:, 1])

        try:

            prot_conc_all.append(
                [float(concentrations[i - 1]) / 1e6 for _ in range(df_temp.shape[0])]
            )  # /1e6 to go from micromolar to molar

        except:

            pass

    mment_factor_all = np.concatenate(mment_factor_all)
    prot_conc_all    = np.concatenate(prot_conc_all)
    signals          = np.concatenate(signals_all).reshape(1, -1)

    return mment_factor_all, signals, prot_conc_all, wavelength_useful