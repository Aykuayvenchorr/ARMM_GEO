class Wellpad:
    def __init__(self, wellpads: [tuple]):
        self.wellpads = wellpads
        self.dict_wellpads : {str: [int, str, bool]} = {}
        self.list_wellpads = []

    def get_dict_wellpads(self):
        # name: [id, luid, rel]
        for wellpad in self.wellpads:
            self.dict_wellpads[wellpad[1]] = [wellpad[0], wellpad[2], wellpad[3]]
        return self.dict_wellpads

    def get_amount_wellpads_for_one_lic_area(self):
        pass

    def create_list_wellpads_for_one_lic_area(self, luid):
        # Добавляем эту строчку здесь, так как в этой функции формируется список имен, которые пойдут в комбобокс
        # self.list_wellpads.append('-------')
        for key, value in self.dict_wellpads.items():
            # добавляем в список площадок для одного ЛУ только актуальные площадки
            # надо, чтобы видно было все площадки, но тут оставлю пока только актуальные
            # if luid in value and value[2]:
            if luid in value:
                self.list_wellpads.append(key)
        return self.list_wellpads











