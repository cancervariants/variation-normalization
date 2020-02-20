class ValidationResult:
    def __init__(self,
            classification,
            is_valid,
            confidence_score,
            concise_description = "",
            human_description = "",
            errors = []):
        self.classification = classification
        self.is_valid = is_valid
        self.confidence_score = confidence_score
        self.concise_description = concise_description
        self.human_description = human_description
        self.errors = errors
