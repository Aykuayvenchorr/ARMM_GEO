class LicArea:
    def __init__(self, areas: [tuple]):
        self.areas = areas
        self.dict_lic_area: {str: str} = {}

    def get_dict_lic_area(self):
        self.dict_lic_area['-------'] = '0000'
        for area in self.areas:
            # Словарь name: id
            self.dict_lic_area[area[1]] = area[0]

        return self.dict_lic_area
