<template>
  <q-page padding>

    <ReactionErrorBanner :store="reactionStore" />

    <div class="row">

      <div class="col-12 col-md-8">
        <Composer />
        <PDFViewer class="q-mt-lg" />
      </div>

      <div class="col-12 col-md-4 q-mt-md">

        <div class="q-pa-md">
          <TextArea
            v-model="reactionStore.textAreaData"
          />
          <Selections :store="reactionStore" />
          <ApplyReset
            :apply="reactionStore.apply"
            :reset="reactionStore.reset"
            :isApplying="reactionStore.isApplying"
            :isResetting="reactionStore.isResetting"
          />
        </div>

      </div>
    </div>

  </q-page>
</template>

<script setup>
import { useReactionStore } from 'src/stores/reaction'
import Composer from 'components/Reaction/Composer.vue'
import PDFViewer from 'components/Reaction/PDFViewer.vue'
import TextArea from 'components/Reaction/TextArea.vue'
import Selections from 'components/Shared/Selections.vue'
import ApplyReset from 'components/Shared/ApplyReset.vue'
import ReactionErrorBanner from 'components/Reaction/ErrorBanner.vue'

import { LocalStorage } from 'quasar'
import { useSettings } from 'src/stores/settings'

const reactionStore = useReactionStore()
const settings = useSettings()
const shouldHidetooltip = settings.filterByTerm("tag", "hideTooltips", "reaction")
reactionStore.hideTooltips(shouldHidetooltip)


</script>
