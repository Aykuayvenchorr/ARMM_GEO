class Wellpad:
    def __init__(self, wellpads: [tuple]):
        self.wellpads = wellpads
        self.dict_wellpads : {str: [int, str, bool]} = {}
        self.list_wellpads = []

    def get_dict_wellpads(self):
        for wellpad in self.wellpads:
            self.dict_wellpads[wellpad[1]] = [wellpad[0], wellpad[2], wellpad[3]]
        return self.dict_wellpads

    def get_amount_wellpads_for_one_lic_area(self):
        pass

    def create_list_wellpads_for_one_lic_area(self, luid):
        for key, value in self.dict_wellpads.items():
            # добавляем в список площадок для одного ЛУ только актуальные площадки
            if luid in value and value[2]:
                self.list_wellpads.append(key)
        return self.list_wellpads


