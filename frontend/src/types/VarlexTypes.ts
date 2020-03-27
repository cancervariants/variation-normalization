export interface TokenResponse {
  searchTerm: string;
  tokens: Token[];
}

export interface Token {
  token: string;
  tokenType: string;
  matchType: string;
  inputString: string;
}

export interface ClassificationResponse {
  searchTerm: string;
  classifications: Classification[];
}

export interface Classification {
  classificationType: string;
  matchingTokens: string[];
  nonMatchingTokens: string[];
  allTokens: Token[];
  confidence: string;
}

export interface ValidationResponse {
  searchTerm: string;
  validationSummary: ValidationSummary;
}

export interface ValidationSummary {
  validResults: ValidationResult[];
  invalidResults: ValidationResult[];
}

export interface ValidationResult {
  classification: Classification;
  isValid: boolean;
  confidenceScore: number;
  humanDescription: string;
  conciseDescription: string;
  errors: string[];
}
