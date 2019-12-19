import axios from "axios";
import { TokenResponse, ClassificationResponse } from "../types/VarlexTypes";
import { varlexApiDomain } from "../services/Config";

const varlexApi = axios.create({
  baseURL: varlexApiDomain(),
  headers: { Accept: "application/json" }
});

export async function getTokens(query: string): Promise<TokenResponse | null> {
  try {
    const resp = await varlexApi.get<TokenResponse>("/tokens/", {
      params: { q: query }
    });

    return resp.data;
  } catch (err) {
    console.log(err);
    return null;
  }
}

export async function getClassifications(
  query: string
): Promise<ClassificationResponse | null> {
  try {
    const resp = await varlexApi.get<ClassificationResponse>(
      "/classifications/",
      {
        params: { q: query }
      }
    );

    return resp.data;
  } catch (err) {
    console.log(err);
    return null;
  }
}
