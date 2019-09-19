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

default_controls = {'minimum combined score': .4, 'minimum neighborhood score': 0, 'minimum gene fusion score': 0, 'minimum cooccurence score': 0, 'minimum coexpression score': 0, 'minimum experimental score': 0, 'minimum database score': 0, 'minimum textmining score': 0, 'minimum best non-textmining component score': 0, 'maximum number of genes': 8}

default_types = {'minimum combined score': 'double', 'minimum neighborhood score': 'double', 'minimum gene fusion score': 'double', 'minimum cooccurence score': 'double', 'minimum coexpression score': 'double', 'minimum experimental score': 'double', 'minimum database score': 'double', 'minimum textmining score': 'double', 'minimum best non-textmining component score': 'double', 'maximum number of genes': 'int'}

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
    limit = get_control(controls, 'maximum number of genes')
    required_score = get_control(controls, 'minimum combined score')
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

        if combined_score >= get_control(controls, 'minimum combined score') and nscore >= get_control(controls, 'minimum neighborhood score') and fscore >= get_control(controls, 'minimum gene fusion score') and pscore >= get_control(controls, 'minimum cooccurence score') and ascore >= get_control(controls, 'minimum coexpression score') and escore >= get_control(controls, 'minimum experimental score') and dscore >= get_control(controls, 'minimum database score') and tscore >= get_control(controls, 'minimum textmining score') and max(nscore,fscore,pscore,ascore,escore,dscore) >= get_control(controls, 'minimum best non-textmining component score'):
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
        name = 'STRING protein-protein interaction',
        function = 'expander',
        #operation = 'interaction',
        #ui_label = 'STRING',
        #source_url = 'https://string-db.org/cgi/help.pl?sessionId=Gp0OeWrM9UF4&subpage=api',
        description = 'Gene-list expander based on STRING protein-protein functional interactions (https://string-db.org/).',
        parameters = [
            Parameter(
                name = 'minimum combined score',
                default = default_controls['minimum combined score'],
                type = default_types['minimum combined score']
            ),
            Parameter(
                name = 'minimum neighborhood score',
                default = default_controls['minimum neighborhood score'],
                type = default_types['minimum neighborhood score']
            ),
            Parameter(
                name = 'minimum gene fusion score',
                default = default_controls['minimum gene fusion score'],
                type = default_types['minimum gene fusion score']
            ),
            Parameter(
                name = 'minimum cooccurence score',
                default = default_controls['minimum cooccurence score'],
                type = default_types['minimum cooccurence score']
            ),
            Parameter(
                name = 'minimum coexpression score',
                default = default_controls['minimum coexpression score'],
                type = default_types['minimum coexpression score']
            ),
            Parameter(
                name = 'minimum experimental score',
                default = default_controls['minimum experimental score'],
                type = default_types['minimum experimental score']
            ),
            Parameter(
                name = 'minimum database score',
                default = default_controls['minimum database score'],
                type = default_types['minimum database score']
            ),
            Parameter(
                name = 'minimum textmining score',
                default = default_controls['minimum textmining score'],
                type = default_types['minimum textmining score']
            ),
            Parameter(
                name = 'minimum best non-textmining component score',
                default = default_controls['minimum best non-textmining component score'],
                type = default_types['minimum best non-textmining component score']
            ),
            Parameter(
                name = 'maximum number of genes',
                default = default_controls['maximum number of genes'],
                type = default_types['maximum number of genes']
            ),
        ],
        required_attributes = ['identifiers.entrez','gene_symbol']
    )
