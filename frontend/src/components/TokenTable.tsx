import React, { Component } from "react";
import { Table, Input, InputOnChangeData } from "semantic-ui-react";
import { TokenResponse } from "../types/VarlexTypes";
import { getTokens } from "../services/VarlexApi"

type State = {
    tokenResponse: TokenResponse | null;
    activeTimeout: number | null;
};

export class TokenTable extends Component<{}, State> {

    state: State = {
        tokenResponse: null,
        activeTimeout: null
    }

    render() {
        return (
            <div>
                <Input icon='search' placeholder='Tokenize' onChange={this.onSearchChanged} />
                <Table celled>
                    <Table.Header>
                        <Table.Row>
                            <Table.HeaderCell>Term</Table.HeaderCell>
                            <Table.HeaderCell>Token Type</Table.HeaderCell>
                        </Table.Row>
                    </Table.Header>
                    {this.tableContents()}
                </Table>
            </div>
        );
    }

    componentDidMount() { this.tokenize() }

    private onSearchChanged = (event: React.ChangeEvent<HTMLInputElement>, data: InputOnChangeData) => {
        if (this.state.activeTimeout) {
            window.clearTimeout(this.state.activeTimeout);
        }

        let searchTerm = event.target.value;
        let newTimer = window.setTimeout(() => { this.tokenize(searchTerm) }, 500);

        this.setState((prevState) => {
            return {
                ...prevState,
                activeTimeout: newTimer
            }
        });
    }

    private tokenize = (searchTerm?: string) => {
        getTokens(searchTerm || '')
            .then(tokenResponse => this.setState({ tokenResponse: tokenResponse }));
    }

    private tableContents = (): JSX.Element => {
        if (this.state.tokenResponse && this.state.tokenResponse.tokens.length > 0) {
            const rows = this.state.tokenResponse.tokens.map(token =>
                <Table.Row key={token.term}>
                    <Table.Cell>{token.term}</Table.Cell>
                    <Table.Cell>{token.type}</Table.Cell>
                </Table.Row>
            );
            return <Table.Body>{rows}</Table.Body>
        } else {
            return <div>No Matches...</div>;
        }
    }
}
