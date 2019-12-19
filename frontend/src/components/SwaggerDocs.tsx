import SwaggerUI from 'swagger-ui-react'
import "swagger-ui-react/swagger-ui.css"
import React, { Component } from 'react'
import { varlexApiDomain } from "../services/Config";


export class SwaggerDocs extends Component<{}, {}> {
    render() {
        let url: string = `${varlexApiDomain()}openapi.json`
        return <SwaggerUI url={url} docExpansion='list' />
    }
}