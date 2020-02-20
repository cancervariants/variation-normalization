class Classification:
    def __init__(self,
        classification_type,
        matching_tokens,
        non_matching_tokens,
        all_tokens,
        confidence
    ):
        self.classification_type = classification_type
        self.matching_tokens = matching_tokens
        self.non_matching_tokens = non_matching_tokens
        self.all_tokens = all_tokens
        self.confidence = confidence

