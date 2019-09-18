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
import requests

species = "9606"

default_controls = {'min combined score': .4, 'min neighborhood score': 0, 'min gene fusion score': 0, 'min phylogenic profile score': 0, 'min coexpression score': 0, 'min experimental score': 0, 'min database score': 0, 'min textmining score': 0, 'min non-textmining score': 0, 'limit': 10}

default_types = {'min combined score': 'double', 'min neighborhood score': 'double', 'min gene fusion score': 'double', 'min phylogenic profile score': 'double', 'min coexpression score': 'double', 'min experimental score': 'double', 'min database score': 'double', 'min textmining score': 'double', 'min non-textmining score': 'double', 'limit': 'int'}

def transform_post(query):  # noqa: E501
    """transform_post

     # noqa: E501

    :param query: Performs transformer query.
    :type query: dict | bytes

    :rtype: List[GeneInfo]
    """
    if connexion.request.is_json:
        query = TransformerQuery.from_dict(connexion.request.get_json())  # noqa: E501

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"

    controls = {control.name:control.value for control in query.controls}
    #Add the originally input genes
    gene_list = query.genes
    genes = dict([(entrez_gene_id(gene) if entrez_gene_id(gene) != None else gene.gene_id, None) for gene in gene_list])
    
    my_genes = genes.keys()
    limit = get_control(controls, 'limit')
    required_score = get_control(controls, 'min combined score')
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
       sys.stderr.print("Error %s\n" % response.status_code)
       return []

    ## Read and parse the results

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

        if combined_score >= get_control(controls, 'min combined score') and nscore >= get_control(controls, 'min neighborhood score') and fscore >= get_control(controls, 'min gene fusion score') and pscore >= get_control(controls, 'min phylogenic profile score') and ascore >= get_control(controls, 'min coexpression score') and escore >= get_control(controls, 'min experimental score') and dscore >= get_control(controls, 'min database score') and tscore >= get_control(controls, 'min textmining score') and max(nscore,fscore,pscore,ascore,escore,dscore) >= get_control(controls, 'min non-textmining score'):
            genes[partner_gene_id] = GeneInfo(
                gene_id = partner_gene_id,
                attributes = [
                   Attribute(
                      name = 'interacting partner',
                      value = query_gene_symbol,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'combined score',
                      value = combined_score,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'nscore',
                      value = nscore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'fscore',
                      value = fscore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'pscore',
                      value = pscore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'ascore',
                      value = ascore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'escore',
                      value = escore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'dscore',
                      value = dscore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'tscore',
                      value = tscore,
                      source = 'STRING interaction expander'
                    ),
                   Attribute(
                      name = 'gene_symbol',
                      value = partner_gene_symbol,
                      source = 'STRING interaction expander'
                    ),
                ]
                #identifiers = GeneInfoIdentifiers(gene_symbol = partner_gene_symbol)
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

def get_control(controls, control):
    value = controls[control] if control in controls else default_controls[control]
    if default_types[control] == 'double':
        return float(value)
    elif default_types[control] == 'Boolean':
        return bool(value)
    elif default_types[control] == 'int':
        return int(value)
    else:
        return value


def transformer_info_get():  # noqa: E501
    """Retrieve transformer info

    Provides information about the transformer. # noqa: E501


    :rtype: TransformerInfo
    """

    return TransformerInfo(
        name = 'STRING interaction expander',
        function = 'expander',
        description = 'Gene-list expander based on STRING interactions.',
        parameters = [
            Parameter(
                name = 'min combined score',
                type = default_controls['min combined score'],
                default = default_types['min combined score']
            ),
            Parameter(
                name = 'min neighborhood score',
                type = default_controls['min neighborhood score'],
                default = default_types['min neighborhood score']
            ),
            Parameter(
                name = 'min gene fusion score',
                type = default_controls['min gene fusion score'],
                default = default_types['min gene fusion score']
            ),
            Parameter(
                name = 'min phylogenic profile score',
                type = default_controls['min phylogenic profile score'],
                default = default_types['min phylogenic profile score']
            ),
            Parameter(
                name = 'min coexpression score',
                type = default_controls['min coexpression score'],
                default = default_types['min coexpression score']
            ),
            Parameter(
                name = 'min experimental score',
                type = default_controls['min experimental score'],
                default = default_types['min experimental score']
            ),
            Parameter(
                name = 'min database score',
                type = default_controls['min database score'],
                default = default_types['min database score']
            ),
            Parameter(
                name = 'min textmining score',
                type = default_controls['min textmining score'],
                default = default_types['min textmining score']
            ),
            Parameter(
                name = 'min non-textmining score',
                type = default_controls['min non-textmining score'],
                default = default_types['min non-textmining score']
            ),
            Parameter(
                name = 'limit',
                type = default_controls['limit'],
                default = default_types['limit']
            ),
        ],
        required_attributes = ['identifiers.entrez','gene_symbol']
    )
