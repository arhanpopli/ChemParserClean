<template>
  <div
    class="text-center q-mx-xl"
    style="width:95%"
  >
    <div class="q-px-xl">
      <q-input
        id="search"
        v-model="molStore.searchTerm"
        dense
        placeholder="Search for a compound by name (e.g aspirin)"
        v-on:keyup.enter="molStore.search"
      >
        <template v-slot:append
            v-if="molStore.isSearching"
          >
            <q-btn
              :loading="molStore.isSearching"
              flat
              round
              color="grey-9"
            />
        </template>
        <template v-slot:append v-else>
          <q-icon
            clickable
            name="fa fa-search"
            size='xs'
            @click="molStore.search"
          />
        </template>
      </q-input>
    </div>

    <div class="q-pa-md">
      <TextArea
        class="full-width"
        rows='25'
        v-model='molStore.textAreaData'
      />
    </div>


    <div class="float-right q-mr-lg">
      <q-btn
        flat
        round
        color="grey-8"
        icon="content_copy"
        @click="molStore.copyToClipBoard"
      >
      </q-btn>
    </div>

    <Tooltip :store="molStore"/>

    <div class="q-mt-lg">
      <SubmitButton
        color="primary"
        label="Submit"
        @click="molStore.submit"
        :loading="molStore.isSubmitting"
      />
    </div>
  </div>
</template>

<script setup>

import { useMoleculeStore } from 'src/stores/molecule'
import SubmitButton from "components/Shared/Button.vue"
import TextArea from 'components/Shared/TextArea.vue'
import Tooltip from 'components/Shared/Tooltip.vue'

const molStore = useMoleculeStore()


</script>
