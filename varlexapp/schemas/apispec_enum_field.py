from marshmallow_enum import EnumField

#Fix for generating OpenAPI enums as documented here:
#https://github.com/justanr/marshmallow_enum/issues/24
#until the pr can be merged
class ApispecEnumField(EnumField):
    def __init__(self, enum, by_value=False, load_by=None, dump_by=None, error='', *args,
            **kwargs):
        super().__init__(enum, by_value, load_by, dump_by, error, *args, **kwargs)
        self.metadata['enum'] = [e.name for e in enum]