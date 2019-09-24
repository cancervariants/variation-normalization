class Token:
    def __init__(self, token, token_type):
        self.token = token
        self.token_type = token_type


    def as_json(self):
        return {
            'term': self.token,
            'type': self.token_type
        }
