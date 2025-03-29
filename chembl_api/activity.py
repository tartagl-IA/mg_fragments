import os

from chembl_webresource_client.settings import Settings

Settings.Instance().FAST_SAVE = False

from typing import Generator

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.query_set import QuerySet

os.environ["CHEMBL_WS_CLIENT_CACHE"] = "0"

client_activity = new_client.activity


def get_activity_from_target_list(
    target_chembl_id_list: list[str],
) -> Generator[QuerySet, None, None]:
    all_results = []
    offset = 0
    batch_size = 1000  # Max recommended limit

    while True:
        batch = client_activity.filter(target_chembl_id__in=target_chembl_id_list)[
            offset : offset + batch_size
        ]

        if not batch:  # Stop when no more data
            break

        # all_results.extend(batch)
        offset += batch_size
        for activity in batch:
            yield activity
