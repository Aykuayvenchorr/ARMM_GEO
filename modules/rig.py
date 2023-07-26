class Rig:
    def __init__(self, rigs: [tuple]):
        self.rigs = rigs
        self.dict_rigs: {str: [int, int, int, float, int, float, bool]} = {}
        self.list_rigs = []
        self.sort_rigs = []

    def get_dict_rigs(self):
        # name: id, type, wellpad_id, nds, nwells, radius, rel, geom
        for rig in self.rigs:
            self.dict_rigs[rig[1]] = [rig[0], rig[2], rig[3], rig[5], rig[6], rig[7], rig[9], rig[8]]
        return self.dict_rigs

    def create_list_rigs_for_one_wellpad(self, wellpad_id):
        # Добавляем эту строчку здесь, так как в этой функции формируется список имен, которые пойдут в комбобокс
        # self.list_rigs.append('-------')
        self.sort_rigs = []
        self.list_rigs = []
        for key, value in self.dict_rigs.items():
            if wellpad_id == value[2]:
                self.list_rigs.append(key)
        self.sort_rigs = sorted(self.list_rigs)
        return self.sort_rigs

    # def null_rigs(self):
    #     self.sort_rigs = []
    #     self.list_rigs = []

