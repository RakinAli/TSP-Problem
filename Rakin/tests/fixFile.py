# Define the input and output file names
input_file_name = "att48.txt"
output_file_name = "att48_updated.txt"

# Open the input file and read the lines
with open(input_file_name, "r") as file:
    lines = file.readlines()

# Open the output file in write mode
with open(output_file_name, "w") as file:
    # Write the number of nodes line as is
    file.write(lines[0])

    # Process each subsequent line
    for line in lines[1:]:
        # Split the line into components
        components = line.split()
        # Write only the x and y positions (the second and third components)
        file.write(components[1] + " " + components[2] + "\n")

print(
    f"The updated file is saved as '{output_file_name}' with only the x and y positions."
)
