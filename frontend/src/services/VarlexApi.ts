import axios from "axios";
import { TokenResponse } from "../types/VarlexTypes";

const varlexApi = axios.create({
  baseURL: "http://localhost:5000",
  headers: { Accept: "application/json" }
});

export async function getTokens(query: string): Promise<TokenResponse | null> {
  try {
    const resp = await varlexApi.get<TokenResponse>("/tokens", {
      params: { q: query }
    });

    return resp.data;
  } catch (err) {
    console.log(err);
    return null;
  }
}
