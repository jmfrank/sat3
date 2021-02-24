

# Class for keeping track of a particular genome and associated files.
class gen:

    def __init__(self, base_dir, genome_fa):
        self.base_dir = base_dir
        self.genome_fa = genome_fa

    # Add
    def add_field(self, field_name, file_str):
        setattr(self, field_name, file_str)
        #self.(field_name) = file_str

