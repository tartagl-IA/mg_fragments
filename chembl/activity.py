from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.query_set import QuerySet

def get_activity_from_target(target_chembl_id) -> QuerySet:
    return new_client.activity.filter(target_chembl_id=target_chembl_id)
