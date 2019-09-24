export interface TokenResponse {
  searchTerm: string;
  tokens: Token[];
}

export interface Token {
  token: string;
  tokenType: string;
}
