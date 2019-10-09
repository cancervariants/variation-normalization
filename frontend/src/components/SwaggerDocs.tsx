import SwaggerUI from 'swagger-ui-react'
import "swagger-ui-react/swagger-ui.css"
import React, { Component } from 'react'


export class SwaggerDocs extends Component<{}, {}> {
    render() {
        return <SwaggerUI url='http://dev.varlex.org/openapi.json' docExpansion='list' />
    }
}