class LicArea:
    def __init__(self, areas: [tuple]):
        self.areas = areas
        self.dict_lic_area: {str: str} = {}

    def get_dict_lic_area(self):
        for area in self.areas:
            self.dict_lic_area[area[1]] = area[0]

        return self.dict_lic_area
