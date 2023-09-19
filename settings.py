import os

plugin_dir = os.path.dirname(__file__)
file_path = os.path.join(
            plugin_dir)

env = {
    'path_doc_save': f'{file_path}/documents/docs',
    'path_target_save': f'{file_path}/documents/targets'
}
