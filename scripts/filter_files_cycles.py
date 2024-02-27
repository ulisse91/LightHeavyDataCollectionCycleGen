import os
import csv
import sys

def filter_csv_file(file_path, substring_to_find, filtered_folder, rename_files):
    filtered_lines = []

    with open(file_path, 'r') as file:
        reader = csv.reader(file)

        try:
            # Copy the first four lines to the filtered list
            for _ in range(4):
                line = next(reader)
                filtered_lines.append(line)

            filtered_count = 0

            # Iterate through the remaining lines
            for row in reader:
                if any(char in row[-1] for char in substring_to_find):
                    filtered_lines.append(row)
                    filtered_count += 1

        except StopIteration:
            print(f"Error: CSV file '{file_path}' does not have the required structure.")
            return None

    # Modify the third line with the total number of lines filtered and the substring_to_find
    filtered_lines[2] = ['K ' + f'{filtered_count}']

    # Create a new file with the filtered lines
    file_dir, file_name = os.path.split(file_path)
    if rename_files:
        filtered_file_name = file_name.replace('.csv', f"_{''.join(substring_to_find)}.csv")
    else:
        filtered_file_name = file_name
    filtered_file_path = os.path.join(filtered_folder, filtered_file_name)

    with open(filtered_file_path, 'w', newline='') as filtered_file:
        writer = csv.writer(filtered_file)
        writer.writerows(filtered_lines)

    return filtered_file_path


# Check if the required number of command-line arguments is provided
if len(sys.argv) < 5:
    print("Error: Please provide the source folder path, destination folder path, substring(s) to find, and rename option as command-line arguments.")
    sys.exit(1)

# Get the source folder path, destination folder path, substring(s) to find, and rename option from the command-line arguments
source_folder_path = sys.argv[1]
destination_folder_path = sys.argv[2]
substring_to_find = sys.argv[3:-1]
rename_files = bool(int(sys.argv[-1]))

# Create the destination folder if it doesn't exist
os.makedirs(destination_folder_path, exist_ok=True)

# Get the list of files to process from the source folder
files_to_process = [filename for filename in os.listdir(source_folder_path) if filename.startswith('cycles_') and filename.endswith('.csv')]

# Calculate the total number of files to be processed
total_files = len(files_to_process)

# Iterate through files in the source folder
for index, filename in enumerate(files_to_process, start=1):
    source_file_path = os.path.join(source_folder_path, filename)

    filtered_file_path = filter_csv_file(source_file_path, substring_to_find, destination_folder_path, rename_files)
    # if filtered_file_path:
    #     print(f"Filtered file created: {filtered_file_path}")

    # Calculate and display progress percentage
    progress = (index / total_files) * 100
    print(f"Progress: {progress:.2f}%\r", end='', flush=True)

print("\nFiltering process completed.")
