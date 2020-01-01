import connexion
import six

from swagger_server.models.gene_info import GeneInfo  # noqa: E501
from swagger_server.models.transformer_info import TransformerInfo  # noqa: E501
from swagger_server.models.transformer_query import TransformerQuery  # noqa: E501
from swagger_server import util
from swagger_server.models.parameter import Parameter
from swagger_server.models.gene_info import GeneInfoIdentifiers
from swagger_server.models.attribute import Attribute

import sys
import json
import requests

class Transformer:

    def __init__(self, variables):
        with open("transformer_info.json",'r') as f:
            self.info = TransformerInfo.from_dict(json.loads(f.read()))
            self.variables = variables
            self.parameters = dict(zip(variables, self.info.parameters))


    def transform(self, query):
        query_controls = {control.name: control.value for control in query.controls}
        controls = {}
        for variable, parameter in self.parameters.items():
            if parameter.name in query_controls:
                controls[variable] = Transformer.get_control(query_controls[parameter.name], parameter)
            else:
                msg = "required parameter '{}' not specified".format(parameter.name)
                return ({ "status": 400, "title": "Bad Request", "detail": msg, "type": "about:blank" }, 400 )

        if self.info.function == 'producer':
            return self.produce(controls)
        if self.info.function == 'expander':
            return self.expand(query.genes, controls)
        if self.info.function == 'filter':
            return self.filter(query.genes, controls)

        return ({ "status": 500, "title": "Internal Server Error", "detail": self.info.function+" not implemented", "type": "about:blank" }, 500 )

    def produce(self, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Producer not implemented", "type": "about:blank" }, 500 )


    def expand(self, query_genes, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Expander not implemented", "type": "about:blank" }, 500 )


    def filter(self, query_genes, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Filter not implemented", "type": "about:blank" }, 500 )


    @staticmethod
    def get_control(value, parameter):
        if parameter.type == 'double':
            return float(value)
        elif parameter.type == 'Boolean':
            return bool(value)
        elif parameter.type == 'int':
            return int(value)
        else:
            return value


species = "9606"

class StringTransformer(Transformer):

    variables = [
        'minimum combined score',
        'minimum neighborhood score',
        'minimum gene fusion score',
        'minimum cooccurence score',
        'minimum coexpression score',
        'minimum experimental score',
        'minimum database score',
        'minimum textmining score',
        'minimum best non-textmining component score',
        'maximum number of genes'
    ]

    def __init__(self):
        super().__init__(self.variables)


    def expand(self, query_genes, controls):

        gene_list = query_genes
        genes = dict([(entrez_gene_id(gene) if entrez_gene_id(gene) != None else gene.gene_id, None) for gene in gene_list])

        string_api_url = "https://string-db.org/api"
        output_format = "tsv-no-header"
        method = "interaction_partners"
        my_genes = genes.keys()

        limit = controls['maximum number of genes']
        required_score = controls['minimum combined score']
        my_app = "sharpener.ncats.io"

        ## Construct the request

        request_url = string_api_url + "/" + output_format + "/" + method + "?"
        request_url += "identifiers=%s" % "%0d".join(my_genes)
        request_url += "&" + "species=" + species
        request_url += "&" + "limit=" + str(limit)
        request_url += "&" + "required_score=" + str(float(required_score) * 1000)
        request_url += "&" + "caller_identity=" + my_app

        response = requests.get(request_url)
        if response.status_code != 200:
           msg = "Call to %s failed (%s)" % (string_api_url, response.status_code)
           return ({ "status": 500, "title": "Internal Server Error", "detail": msg, "type": "about:blank" }, 500 )

        symbol_to_id = {}

        lines = response.text.split('\n')

        for line in lines:
            l = line.strip().split("\t")
            if len(l) < 13:
                continue
            query_gene_symbol = l[2]
            partner_gene_symbol = l[3]
            combined_score = float(l[5])
            nscore = float(l[6])
            fscore = float(l[7])
            pscore = float(l[8])
            ascore = float(l[9])
            escore = float(l[10])
            dscore = float(l[11])
            tscore = float(l[12])

            if query_gene_symbol in symbol_to_id:
                query_gene_id = symbol_to_id[query_gene_symbol]
            else:
                query_gene_id = map_symbol_to_entrez_id(query_gene_symbol)
                symbol_to_id[query_gene_symbol] = query_gene_id

            if partner_gene_symbol in symbol_to_id:
                partner_gene_id = symbol_to_id[partner_gene_symbol]
            else:
                partner_gene_id = map_symbol_to_entrez_id(partner_gene_symbol)
                symbol_to_id[partner_gene_symbol] = partner_gene_id

            if partner_gene_id in genes:
                continue

            if combined_score >= controls['minimum combined score'] and nscore >= controls['minimum neighborhood score'] and fscore >= controls['minimum gene fusion score'] and pscore >= controls['minimum cooccurence score'] and ascore >= controls['minimum coexpression score'] and escore >= controls['minimum experimental score'] and dscore >= controls['minimum database score'] and tscore >= controls['minimum textmining score'] and max(nscore,fscore,pscore,ascore,escore,dscore) >= controls['minimum best non-textmining component score']:
                genes[partner_gene_id] = GeneInfo(
                    gene_id = partner_gene_id,
                    attributes = [
                       Attribute(
                          name = 'interacting partner',
                          value = query_gene_symbol,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'combined score',
                          value = combined_score,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'neighborhood score',
                          value = nscore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'gene fusion score',
                          value = fscore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'cooccurence score',
                          value = pscore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'coexpression score',
                          value = ascore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'experimental score',
                          value = escore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'database score',
                          value = dscore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'textmining score',
                          value = tscore,
                          source = self.info.name
                        ),
                       Attribute(
                          name = 'gene_symbol',
                          value = partner_gene_symbol,
                          source = self.info.name
                        ),
                    ]
                )
                gene_list.append(genes[partner_gene_id])

        return gene_list


def hgnc_gene_id(gene: GeneInfo):
    """
        Return value of the entrez_gene_id attribute
    """
    if (gene.identifiers is not None and gene.identifiers.hgnc is not None):
        if (gene.identifiers.hgnc.startswith('HGNC:')):
            return gene.identifiers.hgnc[5:]
        else:
            return gene.identifiers.hgnc
    return None


def entrez_gene_id(gene: GeneInfo):
    """
        Return value of the entrez_gene_id attribute
    """
    if (gene.identifiers is not None and gene.identifiers.entrez is not None):
        if (gene.identifiers.entrez.startswith('NCBIGene:')):
            return gene.identifiers.entrez[9:]
        else:
            return gene.identifiers.entrez
    return None


def map_symbol_to_entrez_id(gene_symbol):
    request_url = "https://mygene.info/v3/query?q=%s" % gene_symbol
    response = requests.get(request_url)
    if response.status_code != 200:
        return None
    response_object = eval(response.text)
    if 'hits' not in response_object:
        return None
    for hit in response_object['hits']:
        if 'taxid' in hit and str(hit['taxid']) == str(species):
            if 'entrezgene' in hit:
                return "NCBIGene:%s" % (hit['entrezgene'])
    return None



