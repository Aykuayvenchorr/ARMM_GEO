class Rig:
    def __init__(self, rigs: [tuple]):
        self.rigs = rigs
        self.dict_rigs: {int: str} = {}
        self.list_rigs = []

    def get_dict_rigs(self):
        for rig in self.rigs:
            self.dict_rigs[rig[1]] = [rig[3]]
        return self.dict_rigs

    def create_list_rigs_for_one_wellpad(self, wellpad_id):
        for key, value in self.dict_rigs.items():
            if wellpad_id in value:
                self.list_rigs.append(key)
        return self.list_rigs
