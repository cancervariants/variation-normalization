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
  confidence: string;
}
