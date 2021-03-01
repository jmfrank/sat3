

# Class for keeping track of a particular genome and associated files.
class gen:

    def __init__(self, base_dir, genome_fa):
        self.base_dir = base_dir
        self.genome_fa = genome_fa

    # Add
    def add_field(self, field_name, file_str):
        setattr(self, field_name, file_str)

    # Add a file.
    def add_file(self, field_name, file_str,type='relative'):
        if type is 'relative':
            setattr(self,field_name,self.base_dir+file_str)
        elif type is 'full':
            setattr(self, field_name, file_str)


