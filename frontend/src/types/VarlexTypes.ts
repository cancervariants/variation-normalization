export interface TokenResponse {
  searchTerm: string;
  tokens: Token[];
}

export interface Token {
  term: string;
  type: string;
}
